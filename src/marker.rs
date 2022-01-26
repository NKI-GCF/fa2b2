extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate num;
extern crate num_traits;

use crate::kmerconst::KmerConst;
use crate::kmerconst::XmerHash;
use crate::kmerstore::KmerStore;
use crate::new_types::{
    extended_position::{ExtPosEtc, MiniExtPosOri},
    position::{BasePos, Position},
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::to_default::ToDefault;
use crate::xmer_location::XmerLoc;
use anyhow::Result;
use bitvec::prelude::Lsb0;
use crossbeam_channel::Sender;
use noodles_fasta as fasta;

pub struct KmerIter<'a> {
    pub(super) scp: Scope<'a>,
    pub(super) ks: &'a mut KmerStore,
    n_stretch: BasePos,
    goffs: BasePos,
    mini_kmp: Vec<MiniExtPosOri>, //
    pub(crate) repetitive: u32,
    period: Position,
    rep_max_dist: usize,
    tx: Vec<Sender<XmerLoc>>,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub(crate) fn new(ks: &'a mut KmerStore, kc: &'a KmerConst, tx: Vec<Sender<XmerLoc>>) -> Self {
        let scp = Scope::new(kc);
        let rep_max_dist = BasePos::from(ks.rep_max_dist).as_usize();
        let mini_shift = rep_max_dist.next_power_of_two().trailing_zeros();
        KmerIter {
            scp,
            ks,
            n_stretch: BasePos::default(),
            goffs: BasePos::default(),
            mini_kmp: vec![MiniExtPosOri::default(); 1 << mini_shift],
            repetitive: 0,
            period: Position::default(),
            rep_max_dist,
            tx,
        }
    }
    /// clear sequence accounting and scope for the processing of another chromosome
    fn init_chromosome_accounting(&mut self) -> Position {
        // we start with no offset on contig, if starting with N's, the stored goffs gets updated
        self.goffs.to_default();
        self.n_stretch.to_default();
        self.repetitive = 0;
        // TODO: allow repetition to include N-stretch - if both sides of N-stretch show the same repetition.
        // If a repetition ends in an N-stretch, thereafter offset to period
        // may differ or the repetition could be different or entirely gone.
        self.period.to_default();
        self.scp.reset();
        self.scp.get_pos()
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
    /// a kmer is recurrent if it reoccurs within rep_max_dist bases, typically 10_000. This is to
    /// handle e.g. transposons or repetitive DNA, which are not mappable otherwise.
    fn is_recurrent(&self, test: &XmerLoc) -> bool {
        let mut x = 0;
        let repeat_idx = test.idx & (self.mini_kmp.len() - 1);
        let stored = &self.mini_kmp[repeat_idx];
        /* FIXME.
         * if BasePos::from(test.p.pos() - stored.p.pos()).to_usize() < self.rep_max_dist {
            if stored.x() == x {
            } else {
                x += 1;
            }
        }*/
        false
    }
    /*fn update_repetitive(&mut self, pd: Position) {
        // XXX waarom werkt dit op mark??
        self.repetitive += 1;
        let mark = self.scp.get_mark();
        let idx = mark.get_idx();
        let stored = self.ks.kmp[idx];
        if stored.is_set() {
            let mark_pos = mark.p.pos();
            let stored_pos = stored.pos();
            if mark_pos != stored_pos {
                if let Some(dist) = mark_pos.get_if_mark_on_period(stored_pos, pd) {
                    if mark_pos > stored_pos {
                        self.ks.extend_repetitive(mark_pos, dist);
                    } else {
                        dbg_print!(
                            "repetitive occurs before already stored [{:x}] p {} <=> stored {}",
                            idx,
                            mark.p,
                            stored
                        );
                        self.ks.replace_repetitive(stored_pos, mark_pos, dist);
                        // FIXME: moet positie nu niet gezet worden, bij less?
                    }
                }
            }
        } else {
            // XXX this occurs, is it an edge case or a bug?
            dbg_print!("repeat with unset mark.idx (b2 corresponds) ??");
        }
    }*/

    fn spawn_process_coding(&mut self) {
        // TODO: allow repetition to include N-stretch - if both sides of N-stretch show the same repetition.
        // If a repetition ends in an N-stretch, thereafter offset to period
        // may differ or the repetition could be different or entirely gone.
        self.period.to_default();

        self.scp.reset();
        let mut coding_pos = self.scp.get_pos();
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
            // FIXME: outdated usage of ks.b2.
            //let b2 = self.ks.b2_for_pos(coding_pos);

            coding_pos.incr();
        }
        // NB. if we push this before complete_and_update_mark, then the test for repeat
        // is_on_last_contig() does not work. Maybe change this if repeat can be handled
        // while storing sequence.
        self.ks.push_contig(ahead_pos, self.goffs);
    }
    fn is_coding_sequence_pending(&self) -> bool {
        self.ks.get_pos() != self.scp.get_pos()
    }
    fn updated_optimal_xmers_only(&mut self, b: u8) -> Option<XmerLoc> {
        let b2 = match b {
            b'A' | b'a' => bitvec![u8, Lsb0; 0, 0],
            b'C' | b'c' => bitvec![u8, Lsb0; 1, 0],
            b'T' | b't' | b'U' | b'u' => bitvec![u8, Lsb0; 0, 1],
            b'G' | b'g' => bitvec![u8, Lsb0; 1, 1],
            n => {
                if n != b'N' && n != b'n' {
                    dbg_print!("Observed odd base {}, treating as ambiguous..", n);
                }
                // TODO / FIXME rather than excluding when Ns occur, try resolving those instances so
                // the seed selects against these as optima. Iterate over the possible sequences in place of the
                // ambiguous. Could try min..max. account which are the sweetspots. Weight is
                // the nr of positions this resolves.
                //
                // past readlen - kmerlen, or we could overcome
                if !self.n_stretch.is_set() && self.is_coding_sequence_pending() {
                    self.spawn_process_coding();
                }
                self.n_stretch.add_assign(1_u64);
                return None;
            }
        };
        if self.n_stretch.is_set() {
            self.finalize_n_stretch();
        }
        //TODO: make sequence storage optional, to a separate file.
        self.ks.store_b2(&b2);
        if let Some(optimum_xmer) = self
            .scp
            .updated_median_xmer(&b2)
            .filter(|mx| !self.is_recurrent(mx))
        {
            // FIXME: put the other strand in the idx top bits (saves reverse complementing).
            // Also, we need more alternative XmerLoc types for distinction along with traits !!
            return Some(optimum_xmer);
        }
        None
    }
    // TODO: spawn tasks per stretch of non-ambiguous sequence (not N).
    pub(crate) fn markcontig(&mut self, kc: &KmerConst, record: fasta::Record) -> Result<()> {
        let chr_coding_start = self.init_chromosome_accounting();
        self.ks.push_contig(chr_coding_start, self.goffs);
        let ext_bits = usize::try_from(self.tx.len().trailing_zeros())?;
        for b in record.sequence().as_ref().iter() {
            if let Some(mut mark) = self.updated_optimal_xmers_only(*b) {
                mark.idx = kc.hash_and_compress(mark.idx, 0);
                let thread_index = mark.get_thread_index(kc.bitlen, ext_bits);
                self.tx[thread_index].send(mark)?;
            }
        }
        dbg_print!("At end, processing last:");
        if self.n_stretch.is_set() {
            self.finalize_n_stretch();
        } else {
            self.spawn_process_coding();
        }

        if record.name() != "test" {
            if let Some(contig) = self.ks.contig.last() {
                let coding: u64 = BasePos::from(contig.twobit - chr_coding_start).into();
                let n_count = u64::try_from(contig.genomic)?;
                let complex: u64 = coding - self.repetitive as u64;
                let tot: u64 = coding + n_count;
                println!(
                "chromosome {}\tcomplex dna:{} of {}({:.2}%)\trepetitive:{}({:.2}%)\tN-count:{}({:.2}%)\t",
                record.name(),
                complex,
                tot,
                100.0 * complex as f64 / tot as f64,
                self.repetitive,
                100.0 * self.repetitive as f64 / tot as f64,
                n_count,
                100.0 * n_count as f64 / tot as f64,
            );
            }
        }
        Ok(())
    }
    fn dist_if_repetitive(&self, stored_p: ExtPosEtc, min_p: ExtPosEtc) -> Option<Position> {
        // FIXME: the contig check was removed. for the upper bound that makes sense for head
        // scope, as upperbound is pos_max rather than previous
        let stored_pos = stored_p.pos();
        let mark_pos = min_p.pos();
        dbg_assert!(mark_pos > stored_pos);
        if self.ks.is_on_last_contig(stored_pos) {
            let dist = mark_pos - stored_pos;
            if dist < self.ks.rep_max_dist {
                return Some(dist);
            }
        }
        None
    }

    /*fn try_store_mark(&mut self, mark: &mut XmerLoc) -> Result<bool> {
        if self.period.is_set() && self.ks.kmp[mark.idx].is_set() {
            self.ks.kmp[mark.idx].set_repetitive();
            return Ok(false);
        }
        let old_stored_p = self.ks.kmp[mark.idx];

        if old_stored_p.is_zero() {
            self.ks.set_kmp(&mark);
            return Ok(false);
        }
        if old_stored_p.is_replaceable_by(mark.p) {
            if old_stored_p.pos() == mark.p.pos() {
                // set and already mark.p. Leave the bit states.
                return Ok(false);
            }
            self.ks.set_kmp(&mark);
            if old_stored_p.x() == mark.p.x() {
                // same extension means same base k-mer origin. this is a duplicate.
                self.ks.kmp[mark.idx].mark_more_recurs_upseq();
            }
            dbg_print!("{} -> ?", mark);
            mark.p = old_stored_p; // to be extended next
        } else if old_stored_p.extension() == mark.p.extension() {
            // Note: same extension and hash means same k-mer origin: identical k-mer sequence.

            // If a kmer occurs multiple times within an extending readlength (repetition),
            // only the first gets a position. During mapping this should be kept in mind.
            if let Some(dist) = self.dist_if_repetitive(old_stored_p, mark.p) {
                //dbg_assert!(dist < self.p.pos() || self.p.is_zero());
                self.period = dist;
                self.ks.kmp[mark.idx].set_repetitive();
                return Ok(false);
            }
            self.ks.kmp[mark.idx].mark_more_recurs_upseq();
        }
        // collision between hashes, the one in mark.p will be extended and tried again.
        Ok(true)
    }*/
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
        let mut kmi = KmerIter::new(ks, kc, vec![]);
        let definition = fasta::record::Definition::new("test", None);
        let sequence = fasta::record::Sequence::from(seq);
        kmi.markcontig(kc, fasta::Record::new(definition, sequence))
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
