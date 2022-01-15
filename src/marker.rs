extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate num;
extern crate num_traits;

use crate::head_scope::HeadScope;
use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use anyhow::Result;
use noodles_fasta as fasta;

pub struct KmerIter<'a> {
    pub(super) scp: HeadScope<'a>,
    pub(super) ks: &'a mut KmerStore,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub(crate) fn new(ks: &'a mut KmerStore, kc: &'a KmerConst) -> Self {
        let scp = HeadScope::new(kc);
        KmerIter { scp, ks }
    }
    fn init_contig(&mut self) -> Position {
        // we start with no offset on contig, if starting with N's, the stored goffs gets updated
        self.scp.reset_for_new_contig();
        let pos = self.scp.p.pos();
        self.ks.push_contig(pos, BasePos::default());
        pos
    }

    // Dit kan beter, met tasks spawning per stretch of non-ambiguous sequence.
    pub(crate) fn markcontig(&mut self, record: fasta::Record) -> Result<()> {
        let coding_offset = self.init_contig();

        let mut seq = record.sequence().as_ref().iter();

        while let Some(b3) = seq.next().map(|b| self.scp.ascii_to_b3(b)) {
            let pos = self.scp.p.pos();
            if let Some(b2) = b3.as_twobit_if_not_n() {
                // no third bit for A, C, T or G.
                // new sequence is also stored, to enable lookup later.
                if let Some(qb) = self.ks.b2.get_mut(pos.byte_pos()) {
                    *qb |= b2.pos_shift(pos).as_u8();
                }
                self.scp.complete_and_update_mark(b2, self.ks)?;
            } else {
                self.scp.elongate_n_stretch(self.ks, pos);
            }
        }
        self.scp.finalize_n_stretch(self.ks);

        if record.name() != "test" {
            if let Some(contig) = self.ks.contig.last() {
                let coding: u64 = BasePos::from(contig.twobit - coding_offset).into();
                let n_count = u64::try_from(contig.genomic)?;
                let complex: u64 = coding - self.scp.repetitive as u64;
                let tot: u64 = coding + n_count;
                println!(
                "chromosome {}\tcomplex dna:{} of {}({:.2}%)\trepetitive:{}({:.2}%)\tN-count:{}({:.2}%)\t",
                record.name(),
                complex,
                tot,
                100.0 * complex as f64 / tot as f64,
                self.scp.repetitive,
                100.0 * self.scp.repetitive as f64 / tot as f64,
                n_count,
                100.0 * n_count as f64 / tot as f64,
            );
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::KmerConst;
    use crate::new_types::extended_position::ExtPosEtc;
    use crate::new_types::position::BasePos;
    use anyhow::Result;
    const SEQLEN: usize = 250;
    const READLEN: u16 = 6;

    fn process<'a>(ks: &'a mut KmerStore, kc: &'a KmerConst, seq: Vec<u8>) -> Result<()> {
        let mut kmi = KmerIter::new(ks, kc);
        let definition = fasta::record::Definition::new("test", None);
        let sequence = fasta::record::Sequence::from(seq);
        kmi.markcontig(fasta::Record::new(definition, sequence))
    }

    #[test]
    fn test_16n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, &kc, b"NNNNNNNNNNNNNNNN"[..].to_owned())?;

        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, Position::default());
        dbg_assert_eq!(ks.contig[0].genomic, 16.into());
        Ok(())
    }
    #[test]
    fn test_1n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, &kc, b"N"[..].to_owned())?;

        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, Position::default());
        dbg_assert_eq!(ks.contig[0].genomic, 1.into());
        Ok(())
    }
    #[test]
    fn test_1n1c1n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, &kc, b"NCN"[..].to_owned())?;

        dbg_assert_eq!(ks.contig.len(), 2);
        dbg_assert_eq!(ks.contig[0].twobit, Position::default());
        dbg_assert_eq!(ks.contig[0].genomic, 1.into());
        dbg_assert_eq!(ks.contig[1].twobit, Position::from_basepos(1_u64));
        dbg_assert_eq!(ks.contig[1].genomic, 3.into());
        Ok(())
    }
    #[test]
    fn test_17c() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, &kc, b"CCCCCCCCCCCCCCCCC"[..].to_owned())?;
        dbg_assert_eq!(ks.kmp.len(), 128);
        let mut first_pos = ExtPosEtc::from_basepos(kc.kmerlen as u64);
        first_pos.set_repetitive();
        let mut seen = 0;
        for i in 1..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert_eq!(ks.kmp[i], first_pos, "i: {}", i);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_1n18c1n() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned())?;
        }
        let mut first_pos = ExtPosEtc::from_basepos(kc.kmerlen as u64);
        first_pos.set_repetitive();
        let mut seen = 0;
        for i in 0..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert_eq!(ks.kmp[i], first_pos, "i: {}", i);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_1n16c() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        {
            process(&mut ks, &kc, b"NCCCCCCCCCCCCCCCC"[..].to_owned())?;
        }
        let mut seen = 0;
        let mut first_pos = ExtPosEtc::from_basepos(kc.kmerlen as u64);
        first_pos.set_repetitive();
        for i in 0..ks.kmp.len() {
            if ks.kmp[i].is_set() {
                dbg_assert_eq!(ks.kmp[i], first_pos, "i: {}", i);
                seen += 1;
            }
        }
        dbg_assert_eq!(seen, 1);
        Ok(())
    }
    #[test]
    fn test_18at() -> Result<()> {
        let kc = KmerConst::new(SEQLEN, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        {
            process(&mut ks, &kc, b"ATATATATATATATATAT"[..].to_owned())?;
        }
        for i in 0..ks.kmp.len() {
            dbg_assert!(
                ks.kmp[i].is_zero() || !ks.kmp[i].is_dup(),
                "[{}], {:?}",
                i,
                ks.kmp[i]
            );
        }
        Ok(())
    }
}
