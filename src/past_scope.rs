use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{
    extension::Extension,
    position::{PosRange, Position},
    twobit::TwoBit,
    xmer::Xmer,
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::xmer_location::XmerLoc;
use anyhow::{ensure, Result};
use std::fmt;

pub struct PastScope<'a> {
    kc: &'a KmerConst,
    p: ExtPosEtc,
    i: usize,
    mod_i: usize,
    plim: PosRange,
    period: Position,
    mark: XmerLoc,
    d: Vec<Xmer>, // misschien is deze on the fly uit ks te bepalen?
    z: Vec<usize>,
}

impl<'a> PastScope<'a> {
    pub(crate) fn new(kc: &'a KmerConst) -> Self {
        PastScope {
            kc,
            p: ExtPosEtc::default(),
            i: 0,
            mod_i: 0,
            plim: PosRange::default(),
            period: Position::default(),
            mark: XmerLoc::default(),
            d: vec![Xmer::new(); kc.no_kmers],
            z: (0..kc.no_kmers).collect(),
        }
    }
    fn rebuild(&mut self, ks: &KmerStore, p: ExtPosEtc, idx: usize) -> Result<()> {
        let pos = Position::from(p);
        ensure!(pos != Position::default());
        self.plim = ks.get_contig_start_end_for_pos(pos);
        let bound = self.kc.get_kmer_boundaries(pos, self.plim);
        dbg_print!("{}", bound);
        let extension = Extension::from(p);

        self.p = ExtPosEtc::from((extension, bound.lower()));
        self.mark.p = ExtPosEtc::from(extension);

        loop {
            if self.update() {
                // we weten extension op voorhand.
                let i = self.pick_mark();
                if self.d[i].pos != self.mark.p.pos() {
                    let mark = self.d[i].get_hash_and_p(self.kc, self.mark.p.x());
                    self.set_mark(&mark);
                    if p.same_pos_and_ext(self.mark.p) {
                        break;
                    }
                    if self.mark.get_idx() == idx {
                        dbg_print!("idx {:x} observed but for {}, not {}", idx, self.mark.p, p);
                    }
                    // XXX ik zou een assertion hier logischer vinden
                    /*assert!(
                        self.p.pos() < bound.ipper(),
                        "kmer {:x} not observed for {:x} !!",
                        idx,
                        p
                    );*/
                    if self.p.pos() >= bound.upper() {
                        dbg_print!("kmer {:x} not observed for {} !!", idx, p);
                        self.p.clear();
                        break;
                    }
                }
            }
            self.increment(ks.b2_for_pos(self.p.pos(), false));
        }
        Ok(())
    }
    fn is_on_contig(&self, pos: Position) -> bool {
        self.plim.has_in_range(pos)
    }
}

impl<'a> Scope for PastScope<'a> {
    fn dist_if_repetitive(
        &self,
        ks: &KmerStore,
        stored_p: ExtPosEtc,
        min_p: ExtPosEtc,
    ) -> Option<Position> {
        let stored_pos = stored_p.pos();
        let mark_pos = min_p.pos();
        dbg_assert!(mark_pos > stored_pos);
        if self.is_on_contig(stored_pos) {
            let dist = mark_pos - stored_pos;
            if dist < ks.rep_max_dist {
                return Some(dist);
            }
        }
        None
    }
    fn pick_mark(&mut self) -> usize {
        let med = self.kc.no_kmers >> 1;
        let i = self
            .z
            .select_nth_unstable_by(med, |&a, &b| self.d[a].cmp(&self.d[b]))
            .1;
        *i
    }
    /// add twobit to k-mers, update k-mer vec, increment pos and update orientation
    /// true if we have at least one kmer.
    fn update(&mut self) -> bool {
        // XXX: function is hot
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.mod_i];
            self.mod_i += 1;
            if self.mod_i == self.kc.no_kmers {
                self.mod_i = 0;
            }
            self.d[self.mod_i] = old_d;
            // FIXME: why off by one?
            self.d[self.mod_i].pos = self.p.pos();
            true
        } else {
            self.i += 1;
            false
        }
    }
    fn increment(&mut self, b2: TwoBit) {
        // first bit is strand bit, set according to kmer orientation bit.
        self.p.set_ori(self.d[self.mod_i].update(self.kc, b2));
        self.p.incr_pos();
    }
    fn set_mark(&mut self, mark: &XmerLoc) {
        dbg_print!("{} (mark)", mark);
        self.mark = *mark;
    }
}

impl<'a> fmt::Display for PastScope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.unshift_pos() as usize;
        let mp = self.mark.p.unshift_pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat((p - self.kc.read_len) * 5);
        let r = p - mp;
        if r == 0 {
            let x = self.kc.read_len - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.get_idx(), n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.read_len {
            let x = self.kc.read_len - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.get_idx(), n << 2, o, s)
        } else {
            //let l = self.kc.read_len - r - n;
            //let ls = if o {" ".repeat(o) + "|"} else {String::from("")};
            //let rs = if l {String::from("|") + &" ".repeat(l << 2)} else {String::from("")};
            //write!(f, "{2}<{3}|{: ^1$x}|{4}>", self.mark.idx, n << 2, o, ls, rs)
            write!(f, "")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::KmerConst;
    use crate::marker::KmerIter;
    use crate::new_types::position::BasePos;
    use anyhow::Result;
    use noodles_fasta as fasta;
    const SEQLEN: usize = 250;
    const READLEN: u16 = 6;

    #[test]
    fn test_reconstruct1() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        let mut kmi = KmerIter::new(&mut ks, &kc);
        let seq_vec = b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..]
            .to_owned();
        let definition = fasta::record::Definition::new("test", None);
        let sequence = fasta::record::Sequence::from(seq_vec);
        kmi.markcontig(fasta::Record::new(definition, sequence))?;
        let mut scp = PastScope::new(&kc);
        let mut seen = 0;
        for hash in 0..kmi.ks.kmp.len() {
            let p = kmi.ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} p:{:?} ]---", hash, p);
                scp.rebuild(&mut kmi.ks, p, hash)?;
                dbg_assert_eq!(scp.mark.p, p.rep_dup_masked(), "[{}]: {:x}", seen, hash);
                seen += 1;
            }
        }
        dbg_assert_eq!(
            seen,
            42,
            "XXX: the number of seen kmers could change, though"
        );
        Ok(())
    }
    #[test]
    fn test_reconstruct_gs4_all() -> Result<()> {
        // all mappable.
        let seqlen: usize = 8;
        let kc = KmerConst::new(seqlen, 2, 0);
        let mut scp = PastScope::new(&kc);

        for gen in 0..=4_usize.pow(seqlen as u32) {
            let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
            let mut kmi = KmerIter::new(&mut ks, &kc);
            let seq_vec: Vec<_> = (0..seqlen)
                .map(|i| match (gen >> (i << 1)) & 3 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'T',
                    3 => 'G',
                    _ => unreachable!(),
                })
                .collect();
            dbg_print!("-- k: {} rl: {} {:#x} seq:", kc.kmerlen, kc.read_len, gen);
            dbg_print!("{:?}", seq_vec);

            let vv: Vec<u8> = seq_vec.into_iter().map(|c| c as u8).collect();
            let definition = fasta::record::Definition::new("test", None);
            let sequence = fasta::record::Sequence::from(vv);
            kmi.markcontig(fasta::Record::new(definition, sequence))?;
            for hash in 0..kmi.ks.kmp.len() {
                let p = kmi.ks.kmp[hash];
                if p.is_set() {
                    dbg_print!("hash: [{:#x}]: p: {:?}", hash, p);
                    scp.rebuild(&mut kmi.ks, p, hash)?;
                    dbg_assert_eq!(
                        scp.mark.p,
                        p.rep_dup_masked(),
                        "reps: {}",
                        kmi.ks.repeat.len()
                    );
                }
            }
        }
        Ok(())
    }
}
