extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate num;
extern crate num_traits;

use crate::head_scope::HeadScope;
use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::position::{BasePos, Position};
use crate::new_types::twobit::TwoBit;
use crate::rdbg::STAT_DB;
use crate::to_default::ToDefault;
use anyhow::Result;
use noodles_fasta as fasta;

pub struct KmerIter<'a> {
    pub(super) scp: HeadScope<'a>,
    pub(super) ks: &'a mut KmerStore,
    n_stretch: BasePos,
    goffs: BasePos,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub(crate) fn new(ks: &'a mut KmerStore, kc: &'a KmerConst) -> Self {
        let scp = HeadScope::new(kc);
        KmerIter {
            scp,
            ks,
            n_stretch: BasePos::default(),
            goffs: BasePos::default(),
        }
    }
    fn init_chromosome_accounting(&mut self) -> Position {
        // we start with no offset on contig, if starting with N's, the stored goffs gets updated
        // we start with no offset on contig, if starting with N's, the stored goffs gets updated
        self.goffs.to_default();
        self.n_stretch.to_default();
        self.scp.partial_reset_get_pos(true)
    }

    /// Finalize an N-stretch, update stored offsets and prepare to process sequence.
    pub(crate) fn finalize_n_stretch(&mut self) {
        if self.n_stretch.is_set() {
            self.goffs += self.n_stretch;
            dbg_print!("N-stretch of {} (goffs: {})", self.n_stretch, self.goffs);
            self.ks.update_contig_genomic_offset(self.goffs);
            self.n_stretch.to_default();
        }
    }

    fn spawn_process_coding(&mut self) -> Result<()> {
        let mut coding_pos = self.scp.partial_reset_get_pos(false);
        let ahead_pos = self.ks.get_pos();
        self.goffs += BasePos::from(ahead_pos - coding_pos);
        dbg_print!(
            "processing coding: {}-{} (goffs: {})",
            BasePos::from(coding_pos),
            BasePos::from(ahead_pos),
            self.goffs
        );

        // if last spawn is still running, block, otherwise spawn task
        // retrieving twobits from first_pos to pos, and processing them in
        // complete_and_update_mark
        while coding_pos != ahead_pos {
            let b2 = self.ks.b2_for_pos(coding_pos, false);
            self.scp.complete_and_update_mark(self.ks, b2)?;
            coding_pos.incr();
        }
        // NB. if we push this before complete_and_update_mark, then the test for repeat
        // is_on_last_contig() does not work. Maybe change this if repeat can be handled
        // while storing sequence.
        self.ks.push_contig(ahead_pos, self.goffs);
        Ok(())
    }

    // TODO: spawn tasks per stretch of non-ambiguous sequence (not N).
    pub(crate) fn markcontig(&mut self, record: fasta::Record) -> Result<()> {
        let chr_coding_start = self.init_chromosome_accounting();
        self.ks.push_contig(chr_coding_start, self.goffs);

        let mut seq = record.sequence().as_ref().iter();

        while let Some(res) = seq.next().map(TwoBit::from_u8) {
            if let Some(b2) = res {
                // NB. test must occur before store_b2(), which updates position.
                if !self.scp.is_coding_sequence_pending(&self.ks) {
                    self.finalize_n_stretch();
                }
                // TODO: also skip repetitive !! then administration should go here!
                self.ks.store_b2(b2)?;
            } else {
                // past readlen - kmerlen, or we could overcome
                if self.scp.is_coding_sequence_pending(&self.ks) {
                    self.spawn_process_coding()?;
                }
                self.n_stretch.add_assign(1_u64);
            }
        }
        dbg_print!("At end, processing last:");
        if self.scp.is_coding_sequence_pending(&self.ks) {
            self.spawn_process_coding()?;
        } else {
            self.finalize_n_stretch();
        }

        if record.name() != "test" {
            if let Some(contig) = self.ks.contig.last() {
                let coding: u64 = BasePos::from(contig.twobit - chr_coding_start).into();
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
        dbg_assert_eq!(ks.contig[1].twobit, Position::from_basepos(1_u64)); //XXX
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
