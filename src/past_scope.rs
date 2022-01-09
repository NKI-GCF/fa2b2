use crate::kmer::Kmer;
use crate::kmer::TwoBit;
use crate::kmerconst::KmerConst;
use crate::kmerloc::{ExtPosEtc, KmerLoc};
use crate::kmerstore::KmerStore;
use crate::new_types::extension::Extension;
use crate::new_types::position::Position;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
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
    d: Vec<Kmer<u64>>, // misschien is deze on the fly uit ks te bepalen?
    z: Vec<usize>,
}

impl<'a> PastScope<'a> {
    pub fn new(ks: &KmerStore, kc: &'a KmerConst, p: ExtPosEtc, idx: usize) -> Result<Self> {
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
            d: vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers],
            z: (0..kc.no_kmers).collect(),
        };
        let x = scp.p.x();
        let bin = kc.get_kmers(x);

        loop {
            if scp.increment(ks.b2_for_p(scp.p.pos(), false).unwrap()) {
                // we weten extension op voorhand.
                let base = scp.get_i() - kc.kmerlen;
                if (scp.set_if_optimum(x, base, bin) && scp.all_kmers()) || scp.remark(false)? {
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
    fn is_before_end_of_contig(&self) -> bool {
        self.p.basepos_to_pos() < self.plim.1
    }
    fn is_on_contig(&self, pos: Position) -> bool {
        pos >= self.plim.0 && pos < self.plim.1
    }
}

impl<'a> Scope for PastScope<'a> {
    fn get_mark(&self) -> Option<&KmerLoc> {
        if self.mark.is_set() {
            Some(&self.mark)
        } else {
            None
        }
    }
    fn get_kc(&self) -> &KmerConst {
        self.kc
    }
    fn get_p(&self) -> ExtPosEtc {
        self.p
    }
    fn get_i(&self) -> usize {
        self.i
    }
    fn get_d(&self, i: usize) -> &Kmer<u64> {
        &self.d[i]
    }
    fn is_repetitive(&self) -> bool {
        self.period.is_set()
    }
    fn clear_p_extension(&mut self) {
        self.p.clear_extension();
    }

    fn increment_for_extension(&mut self, ks: &KmerStore) -> Result<()> {
        let b2 = ks.b2_for_p(self.get_p().pos(), false)?;
        ensure!(
            self.is_before_end_of_contig(),
            "increment_for_extension() runs into end of contig"
        );
        dbg_assert!(self.increment(b2));
        Ok(())
    }

    fn dist_if_repetitive(
        &self,
        stored_p: ExtPosEtc,
        mark_p: ExtPosEtc,
        max_dist: Position,
    ) -> Option<Position> {
        let stored_pos = stored_p.pos();
        let mark_pos = mark_p.pos();
        dbg_assert!(mark_pos > stored_pos);
        if self.is_on_contig(stored_pos) {
            let dist = mark_pos - stored_pos;
            if dist < max_dist {
                return Some(dist);
            }
        }
        None
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
        }
        // first bit is strand bit, set according to kmer orientation bit.
        self.p.set_ori(self.d[self.mod_i].update(b2));
        self.p.incr_pos();
        self.i += 1;
        self.i >= self.kc.kmerlen
    }
    fn extend_p(&mut self) {
        self.p.extend();
    }
    fn mark_reset(&mut self) {
        self.mark.reset();
    }
    fn set_mark(&mut self, idx: usize, p: ExtPosEtc, x: usize) {
        dbg_print!("[{:x}] = {:?} | x({})", idx, p, x);
        self.mark.set(idx, p, x);
    }
    fn set_period(&mut self, period: Position) {
        dbg_assert!(period < self.p.pos());
        self.period = period;
    }
    fn unset_period(&mut self) {
        self.set_period(Position::zero());
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
        kmi.ks.pos_max = (seq_vec.len() as u64).basepos_to_pos();
        let definition = fasta::record::Definition::new("test", None);
        let sequence = fasta::record::Sequence::from(seq_vec);
        kmi.markcontig::<u64>(fasta::Record::new(definition, sequence))?;
        let mut seen = 0;
        for hash in 0..kmi.ks.kmp.len() {
            let p = kmi.ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} p:{:x} ]---", hash, p);
                let scp = PastScope::new(&kmi.ks, &kc, &p, hash)?;
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
            kmi.ks.pos_max = (seqlen as u64).basepos_to_pos();
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
            kmi.markcontig::<u64>(fasta::Record::new(definition, sequence))?;
            for hash in 0..kmi.ks.kmp.len() {
                let p = kmi.ks.kmp[hash];
                if p.is_set() {
                    dbg_print!("hash: [{:#x}]: p: {:#x}", hash, p);
                    let scp = PastScope::new(&kmi.ks, &kc, &p, hash)?;
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
