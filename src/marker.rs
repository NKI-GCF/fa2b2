extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate num;
extern crate num_traits;

use crate::kmerconst::KmerConst;
use std::{cmp, slice::Iter};

use crate::kmerloc::PriExtPosOri;
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::Result;

pub struct KmerIter<'a> {
    n_stretch: u64,
    goffs: u64,
    pub(super) scp: Scope<'a>,
    pub(super) ks: &'a mut KmerStore<u64>,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore<u64>, kc: &'a KmerConst) -> Self {
        let scp = Scope::new((0, u64::max_value()), kc, 0);
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
        self.scp.p.clear_extension();
        self.scp.mark.reset();
        self.scp.plim.0 = self.scp.p.pos();
        self.scp.mod_i = 0;
        // If a repetition ends in an N-stretch, thereafter offset to period
        // may differ or the repetition could be different or entirely gone.
        // TODO: allow repetition to include N-stretch - if both sides of N-stretch show the same repetition.
        self.scp.period = 0;
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
                self.ks.p_max = p.pos() + 4;
                if self.n_stretch > 0 {
                    self.finalize_n_stretch();
                } else if self.scp.period != 0 && dbg_dump_if!(self.scp.mark.is_set(), false) {
                    // XXX self.scp.mark.is_set() can be false here, it seems.
                    let pd = self.scp.period;
                    if self.ks.b2_for_p(self.scp.p - pd)? == b2 {
                        let idx = self.scp.mark.get_idx();
                        let stored = self.ks.kmp[idx];
                        let dist = match self.scp.mark.p.pos().cmp(&stored.pos()) {
                            cmp::Ordering::Greater => self.scp.mark.p.pos() - stored.pos(),
                            cmp::Ordering::Less => {
                                dbg_print!("repetitive occurs before non-repetitive?");
                                stored.pos() - self.scp.mark.p.pos()
                            }
                            cmp::Ordering::Equal => {
                                // FIXME: waarom gebeurt dit? TODO: ignore het niet.
                                //dbg_panic!("revisit [{:x}] {:x} (pd: {})?", idx, stored, pd);
                                dbg_print!("revisit [{:x}] {:x} (pd: {})?", idx, stored, pd);
                                dbg_dump!();
                                pd - 1
                            }
                        };
                        if dist % pd == 0 {
                            repetitive += 1;
                            self.ks.extend_repetitive(idx, dist as u32);
                        }
                    } else {
                        self.scp.period = 0;
                    }
                }
                // scp funcs also used for scope rebuild, therefore ext is set here.
                self.scp.complete_and_update_mark::<u64>(b2, &mut self.ks)?;
            } else {
                if self.scp.i != 0 {
                    dbg_print!("started N-stretch at {}.", p);
                    self.goffs += self.scp.i as u64;
                    self.ks.push_contig(p, self.goffs);

                    // clear all except orientation and position to rebuild at the start of a new contig.
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
        let first_pos = 1 | (kc.kmerlen as u64) << 1;
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
        let first_pos = 1 | (kc.kmerlen as u64) << 1;
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
        let first_pos = 1 | (kc.kmerlen as u64) << 1;
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
    #[test]
    fn test_reconstruct1() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut ks = KmerStore::<u64>::new(kc.bitlen);
        let mut kmi = KmerIter::new(&mut ks, &kc);
        let seq_vec = b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..]
            .to_owned();
        let mut seq = seq_vec.iter();
        kmi.markcontig::<u64>("test", &mut seq)?;
        let mut seen = 0;
        for hash in 0..kmi.ks.kmp.len() {
            let p = kmi.ks.kmp[hash];
            if p.is_set() {
                dbg_print!("---[ {:#x} ]---", hash);
                let scp = Scope::rebuild(&kmi.ks, &kc, &p, hash)?;
                dbg_assert_eq!(scp.mark.p, p, "[{}]: {:x}", seen, hash);
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
        let kc = KmerConst::new(seqlen, 1000);

        for gen in 0..=4_usize.pow(seqlen as u32) {
            let mut ks = KmerStore::<u64>::new(kc.bitlen);
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
            dbg_print!("-- k: {} rl: {} {:#x} seq:", kc.kmerlen, kc.venster, gen);
            dbg_print!("{:?}", seq_vec);

            let vv: Vec<u8> = seq_vec.into_iter().map(|c| c as u8).collect();
            kmi.markcontig::<u64>("test", &mut vv.iter())?;
            for hash in 0..kmi.ks.kmp.len() {
                let p = kmi.ks.kmp[hash];
                if p.is_set() {
                    dbg_print!("hash: [{:#x}]: p: {:#x}", hash, p);
                    let scp = Scope::rebuild(&kmi.ks, &kc, &p, hash)?;
                    dbg_assert_eq!(scp.mark.p, p, "reps: {}", kmi.ks.repeat.len());
                }
            }
        }
        Ok(())
    }
}
