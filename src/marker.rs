extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate num;
extern crate num_traits;

use crate::kmerconst::KmerConst;
use std::{cmp, slice::Iter};

use crate::head_scope::HeadScope;
use crate::kmerloc::PriExtPosOri;
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use anyhow::Result;

pub struct KmerIter<'a> {
    n_stretch: u64,
    goffs: u64,
    pub(super) scp: HeadScope<'a>,
    pub(super) ks: &'a mut KmerStore<u64>,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst) -> Self {
        let scp = HeadScope::new((0, u64::max_value()), kc, 0);
        KmerIter {
            n_stretch: 0,
            goffs: 0,
            scp,
            ks,
        }
    }

    /// sluit een n-stretch af, die we aan het verlengen waren.
    fn finalize_n_stretch(&mut self) {
        dbg_print!("added new contig. Ns:{}", self.n_stretch);
        self.ks.offset_contig(self.n_stretch);
        self.goffs += self.n_stretch;
        self.n_stretch = 0;
        // clear all except orientation and position to rebuild at the start of a new contig.
        if self.scp.p.is_set() {
            // clearing extension on first position fails assert on is_set()
            self.scp.p.clear_extension();
        }
        self.scp.mark.reset();
        self.scp.plim.0 = self.scp.p.pos();
        self.scp.mod_i = 0;
        // If a repetition ends in an N-stretch, thereafter offset to period
        // may differ or the repetition could be different or entirely gone.
        // TODO: allow repetition to include N-stretch - if both sides of N-stretch show the same repetition.
        self.scp.period = 0;
    }

    fn extend_repetive_on_cycle(&mut self, mark_pos: u64, dist: u64, pd: u64) {
        if dist % pd == 0 {
            self.ks.extend_repetitive(mark_pos, dist as u32);
        }
    }

    fn replace_repetive_on_cycle(&mut self, old_pos: u64, new_pos: u64, dist: u64, pd: u64) {
        if dist % pd == 0 {
            self.ks.replace_repetitive(old_pos, new_pos, dist as u32);
        }
    }

    pub fn markcontig<T: PriExtPosOri>(&mut self, chrname: &str, seq: &mut Iter<u8>) -> Result<()> {
        self.goffs = 0;
        self.ks.push_contig(self.scp.p.pos(), self.goffs);

        let mut repetitive = 0_u64;
        let mut n_count = 0_u64;
        let mut tot = 0_u64;

        while let Some(b2) = seq.next().map(|&c| {
            let b2 = (c >> 1) & 0x7;
            // dbg_print!("[{}, {}]: {:x}", self.scp.p.pos() >> 1, c as char, b2);
            dbg_print!("{}: {:x}", c as char, b2);
            b2
        }) {
            tot += 1;
            let p = self.scp.p;
            if b2 < 4 {
                // new sequence is also stored, to enable lookup later.
                if let Some(qb) = self.ks.b2.get_mut(p.byte_pos()) {
                    *qb |= b2 << (p & 6);
                }
                self.ks.p_max = p.pos() + 2;
                if self.n_stretch > 0 {
                    self.finalize_n_stretch();
                } else if self.scp.period != 0 && dbg_dump_if!(self.scp.mark.is_set(), false) {
                    // XXX self.scp.mark.is_set() can be false here, it seems.
                    let pd = self.scp.period;
                    if self.ks.b2_for_p(self.scp.p - pd, Some("(R)"))? == b2 {
                        let idx = self.scp.mark.get_idx();
                        let stored = self.ks.kmp[idx];
                        repetitive += 1;
                        if stored.is_set() {
                            let mark_pos = self.scp.mark.p.pos();
                            let stored_pos = stored.pos();
                            match mark_pos.cmp(&stored_pos) {
                                cmp::Ordering::Greater => {
                                    let dist = mark_pos - stored_pos;
                                    self.extend_repetive_on_cycle(mark_pos, dist, pd);
                                }
                                cmp::Ordering::Less => {
                                    dbg_print!("repetitive occurs before already stored [{:x}] p {:#x} <=> stored {:#x}", idx, self.scp.mark.p, stored);
                                    let _x = self.scp.p.x();
                                    dbg_assert!(_x > 0);
                                    let dist = stored_pos - mark_pos;
                                    self.replace_repetive_on_cycle(stored_pos, mark_pos, dist, pd);
                                    // FIXME: moet positie nu niet gezet worden, bij less?
                                }
                                cmp::Ordering::Equal => {
                                    dbg_print!("minimum remained [{:x}] {:#x}", idx, stored)
                                }
                            }
                        } else {
                            // XXX this occurs, is it an edge case or a bug?
                            dbg_print!("repeat with unset mark.idx (b2 corresponds) ??");
                        }
                    } else {
                        self.scp.period = 0;
                    }
                }
                self.scp.complete_and_update_mark::<u64>(b2, self.ks)?;
            } else {
                if self.scp.i != 0 {
                    dbg_print!("started N-stretch at {}.", p);
                    self.goffs += self.scp.i as u64;
                    self.ks.push_contig(p, self.goffs);

                    self.n_stretch = 0;
                    self.scp.i = 0;
                }
                self.n_stretch += 1;
                n_count += 1;
            }
        }
        if self.n_stretch > 0 {
            dbgx!(self.finalize_n_stretch());
        }
        if chrname != "test" {
            let complex = tot - repetitive - n_count;
            println!(
                "chromosome {}\tcomplex dna:{} of {}({:.2}%)\trepetitive:{}({:.2}%)\tN-count:{}({:.2}%)\t",
                chrname,
                complex,
                tot,
                100.0 * (tot - repetitive - n_count) as f64 / tot as f64,
                repetitive,
                100.0 * repetitive as f64 / tot as f64,
                n_count,
                100.0 * n_count as f64 / tot as f64,
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::KmerConst;
    use anyhow::Result;
    const SEQLEN: usize = 250;

    fn process<'a>(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst, seq: Vec<u8>) -> Result<()> {
        let mut kmi = KmerIter::new(ks, kc);
        kmi.markcontig::<u64>("test", &mut seq.iter())
    }

    #[test]
    fn test_16n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NNNNNNNNNNNNNNNN"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, 0);
        dbg_assert_eq!(ks.contig[0].genomic, 16);
        Ok(())
    }
    #[test]
    fn test_1n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"N"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, 0);
        dbg_assert_eq!(ks.contig[0].genomic, 1);
        Ok(())
    }
    #[test]
    fn test_1n1c1n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCN"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.contig.len(), 2);
        dbg_assert_eq!(ks.contig[0].twobit, 0);
        dbg_assert_eq!(ks.contig[0].genomic, 1);
        dbg_assert_eq!(ks.contig[1].twobit, 2);
        dbg_assert_eq!(ks.contig[1].genomic, 3);
        Ok(())
    }
    #[test]
    fn test_17c() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"CCCCCCCCCCCCCCCCC"[..].to_owned())?;
        }
        dbg_assert_eq!(ks.kmp.len(), 128);
        let mut first_pos = 1 | (kc.kmerlen as u64) << 1;
        first_pos.set_repetitive();
        let mut seen = 0;
        for i in 1..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert!(ks.kmp[i] == first_pos, "[{:x}], {:x}", i, ks.kmp[i]);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_1n18c1n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned())?;
        }
        let mut first_pos = 1 | (kc.kmerlen as u64) << 1;
        first_pos.set_repetitive();
        let mut seen = 0;
        for i in 1..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert!(ks.kmp[i] == first_pos, "[{:x}], {:x}", i, ks.kmp[i]);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_1n16c() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCC"[..].to_owned())?;
        }
        let mut seen = 0;
        let mut first_pos = 1 | (kc.kmerlen as u64) << 1;
        first_pos.set_repetitive();
        for i in 1..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert!(ks.kmp[i] == first_pos, "[{:x}], {:x}", i, ks.kmp[i]);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_18at() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        {
            process(&mut ks, &kc, b"ATATATATATATATATAT"[..].to_owned())?;
        }
        for i in 0..ks.kmp.len() {
            dbg_assert!(
                ks.kmp[i].is_no_pos() || !ks.kmp[i].is_dup(),
                "[{}], {:x}",
                i,
                ks.kmp[i]
            );
        }
        Ok(())
    }
}
