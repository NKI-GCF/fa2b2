use crate::kmerconst::KmerConst;
use crate::kmerloc::{ExtPosEtc, KmerLoc};
use crate::kmerstore::KmerStore;
use crate::new_types::{extension::Extension, position::Position, twobit::TwoBit, xmer::Xmer};
use crate::rdbg::STAT_DB;
use crate::scope::{Scope, WritingScope};
use anyhow::{ensure, Result};
use std::fmt;

pub struct PastScope<'a> {
    kc: &'a KmerConst,
    p: ExtPosEtc,
    i: usize,
    mod_i: usize,
    plim: (Position, Position),
    period: Position,
    mark: KmerLoc,
    d: Vec<Xmer>, // misschien is deze on the fly uit ks te bepalen?
    z: Vec<usize>,
}

impl<'a> PastScope<'a> {
    pub fn new(ks: &mut KmerStore, kc: &'a KmerConst, p: ExtPosEtc, idx: usize) -> Result<Self> {
        let pos = Position::from(p);
        ensure!(pos != Position::zero());
        let plim = ks.get_contig_start_end_for_p(pos);
        let bound = kc.get_kmer_boundaries(pos, plim);
        dbg_print!("{:?}", bound);
        let extension = Extension::from(p);

        let mut scp = PastScope {
            kc,
            p: ExtPosEtc::from((extension, bound.0)),
            i: 0,
            mod_i: 0,
            plim,
            period: Position::zero(),
            mark: KmerLoc::new(usize::max_value(), ExtPosEtc::from(extension)),
            d: vec![Xmer::new(kc.kmerlen as u32); kc.no_kmers],
            z: (0..kc.no_kmers).collect(),
        };

        loop {
            if scp.increment(ks.b2_for_p(scp.p.pos(), false).unwrap()) {
                // we weten extension op voorhand.
                //
                // we weten extension op voorhand.
                let i = scp.pick_mark();
                if scp.d[i].pos != scp.mark.p.pos() {
                    scp.store_mark(ks, i)?;
                    if p.same_pos_and_ext(scp.mark.p) {
                        break;
                    }
                    if scp.mark.get_idx() == idx {
                        dbg_print!(
                            "idx {:x} observed but for {:?}, not {:?}",
                            idx,
                            scp.mark.p,
                            p
                        );
                    }
                    // XXX ik zou een assertion hier logischer vinden
                    /*assert!(
                        scp.p.pos() < bound.1,
                        "kmer {:x} not observed for {:x} !!",
                        idx,
                        p
                    );*/
                    if scp.p.pos() >= bound.1 {
                        dbg_print!("kmer {:x} not observed for {:?} !!", idx, p);
                        scp.p.clear();
                        break;
                    }
                }
            }
        }
        Ok(scp)
    }
    fn is_on_contig(&self, pos: Position) -> bool {
        pos >= self.plim.0 && pos < self.plim.1
    }
}

impl<'a> WritingScope for PastScope<'a> {
    fn set_period(&mut self, period: Position) {
        dbg_assert!(period < self.p.pos());
        self.period = period;
    }
    fn unset_period(&mut self) {
        self.set_period(Position::zero());
    }
    fn is_repetitive(&self) -> bool {
        self.period.is_set()
    }
}

impl<'a> Scope for PastScope<'a> {
    fn get_kc(&self) -> &KmerConst {
        self.kc
    }
    fn get_d(&self, i: usize) -> &Xmer {
        &self.d[i]
    }

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
    fn increment(&mut self, b2: TwoBit) -> bool {
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
        }
        // first bit is strand bit, set according to kmer orientation bit.
        self.p.set_ori(self.d[self.mod_i].update(b2));
        self.p.incr_pos();
        self.i += 1;
        self.i >= self.kc.kmerlen
    }
    fn set_mark(&mut self, idx: usize, p: ExtPosEtc) {
        dbg_print!("[{:x}] = {:?}", idx, p);
        self.mark.set(idx, p);
    }
}

impl<'a> fmt::Display for PastScope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.unshift_pos() as usize;
        let mp = self.mark.p.unshift_pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat((p - self.kc.venster) * 5);
        let r = p - mp;
        if r == 0 {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.get_idx(), n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.venster {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.get_idx(), n << 2, o, s)
        } else {
            //let l = self.kc.venster - r - n;
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

    #[test]
    fn test_reconstruct1() -> Result<()> {
        let kc = KmerConst::new(SEQLEN);
        let mut ks = KmerStore::new(kc.bitlen, 10_000);
        let mut kmi = KmerIter::new(&mut ks, &kc);
        let seq_vec = b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..]
            .to_owned();
        kmi.ks.pos_max = Position::from(BasePos::from(seq_vec.len()));
        let definition = fasta::record::Definition::new("test", None);
        let sequence = fasta::record::Sequence::from(seq_vec);
        kmi.markcontig(fasta::Record::new(definition, sequence))?;
        let mut seen = 0;
        for hash in 0..kmi.ks.kmp.len() {
            let p = kmi.ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} p:{:?} ]---", hash, p);
                let scp = PastScope::new(&mut kmi.ks, &kc, p, hash)?;
                dbg_assert_eq!(scp.mark.p, p.rep_dup_masked(), "[{}]: {:x}", seen, hash);
                seen += 1;
            }
        }
        dbg_assert_eq!(
            seen,
            23,
            "XXX: the number of seen kmers could change, though"
        );
        Ok(())
    }
    #[test]
    fn test_reconstruct_gs4_all() -> Result<()> {
        // all mappable.
        let seqlen: usize = 8;
        let kc = KmerConst::new(seqlen);

        for gen in 0..=4_usize.pow(seqlen as u32) {
            let mut ks = KmerStore::new(kc.bitlen, 10_000);
            let mut kmi = KmerIter::new(&mut ks, &kc);
            kmi.ks.pos_max = Position::from(BasePos::from(seqlen));
            let seq_vec: Vec<_> = (0..seqlen)
                .map(|i| match (gen >> (i << 1)) & 3 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'T',
                    3 => 'G',
                    _ => unreachable!(),
                })
                .collect();
            dbg_print!("-- k: {} rl: {} {:#x} seq:", kc.kmerlen, kc.venster, gen);
            dbg_print!("{:?}", seq_vec);

            let vv: Vec<u8> = seq_vec.into_iter().map(|c| c as u8).collect();
            let definition = fasta::record::Definition::new("test", None);
            let sequence = fasta::record::Sequence::from(vv);
            kmi.markcontig(fasta::Record::new(definition, sequence))?;
            for hash in 0..kmi.ks.kmp.len() {
                let p = kmi.ks.kmp[hash];
                if p.is_set() {
                    dbg_print!("hash: [{:#x}]: p: {:?}", hash, p);
                    let scp = PastScope::new(&mut kmi.ks, &kc, p, hash)?;
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
