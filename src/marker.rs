extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate num;
extern crate num_traits;

use crate::kmerconst::KmerConst;
use crate::kmerconst::XmerHash;
use crate::kmerstore::KmerStore;
use crate::new_types::{
    extended_position::ExtPosEtc,
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
    mini_kmp: Vec<ExtPosEtc>, //
    pub(crate) repetitive: u32,
    period: Position,
    tx: Vec<Sender<XmerLoc>>,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub(crate) fn new(
        ks: &'a mut KmerStore,
        kc: &'a KmerConst,
        tx: Vec<Sender<XmerLoc>>,
    ) -> Result<Self> {
        let scp = Scope::new(kc)?;
        let mini_shift = BasePos::from(ks.rep_max_dist)
            .as_usize()
            .next_power_of_two()
            .trailing_zeros();
        Ok(KmerIter {
            scp,
            ks,
            n_stretch: BasePos::default(),
            goffs: BasePos::default(),
            mini_kmp: vec![ExtPosEtc::default(); 1 << mini_shift],
            repetitive: 0,
            period: Position::default(),
            tx,
        })
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
    fn finalize_n_stretch(&mut self) {
        self.goffs += self.n_stretch;
        dbg_print!("N-stretch of {} (goffs: {})", self.n_stretch, self.goffs);
        self.ks.update_contig_genomic_offset(self.goffs);
        self.n_stretch.to_default();
    }

    // FIXME: the repetitive bit is no longer set, if we just filter repetitive on period.
    // We should but how? send mark with samepose and repetitive bit set? set these only after completion?
    // FIXME, currently this filters every recurrent, but we should only filter if there is at
    // least a completed cycle of repeats.

    /// a kmer is recurrent if it reoccurs within rep_max_dist bases, typically 10_000. This is to
    /// handle e.g. transposons or repetitive DNA, which are not mappable otherwise.
    fn is_recurrent(&mut self, mut median_xmer: XmerLoc) -> Option<XmerLoc> {
        let mask = self.mini_kmp.len() - 1;
        let mut mx_clone = median_xmer;
        let mut ret = Some(median_xmer);
        let mut test = &mut median_xmer;
        loop {
            let repeat_idx = test.idx & mask;
            let stored = self.mini_kmp[repeat_idx];

            if stored.is_zero() || test.p.pos() > self.ks.rep_max_dist + stored.pos() {
                dbg_assert!(
                    stored.pos() + Position::from_basepos((self.scp.kc.no_kmers / 2) as u64)
                        < test.p.pos(),
                    "{} >= {} ?",
                    stored.pos() + Position::from_basepos((self.scp.kc.no_kmers / 2) as u64),
                    test.p.pos()
                );
                self.mini_kmp[repeat_idx].set(test.p);
                // unsetting dup isn't necessary, if from a stored, we're a dup anyway.
                break;
            }
            // TODO: allow one long stretch of Ns in between?
            if self.ks.is_on_last_contig(stored.pos()) {
                if stored.gt_ext_or_eq_and_ge_pos(test.p) {
                    if stored.is_dup() {
                        ret = None;
                    }
                    if stored.pos() == test.p.pos() {
                        break;
                    }
                    self.mini_kmp[repeat_idx].set(test.p);
                    if stored.x() == test.p.x() {
                        // same extension means same base k-mer origin. this is a duplicate.
                        self.mini_kmp[repeat_idx].set_dup();
                    }
                    test = &mut mx_clone;
                    test.p.set(stored); // to be extended next
                } else if stored.extension() == test.p.extension() {
                    // Note: same extension and hash means same k-mer origin: identical k-mer sequence.
                    if stored.is_dup() {
                        ret = None;
                    }
                    self.mini_kmp[repeat_idx].set_dup();
                }
                // From the XmerHash trait, extension fails on last, when all ext bits are set.
                if self.extend_xmer(test).is_err() {
                    dbg_print!("couldn't extend {}", test);
                    break;
                }
            }
        }
        ret
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

    fn finalize_coding(&mut self) {
        // TODO: allow repetition to include N-stretch - if both sides of N-stretch show the same repetition.
        // If a repetition ends in an N-stretch, thereafter offset to period
        // may differ or the repetition could be different or entirely gone.
        self.period.to_default();

        self.scp.reset();

        self.ks.push_contig(self.ks.get_pos(), self.goffs);
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
                // the seed selects against these as optimal. Iterate over the possible sequences in place of the
                // ambiguous. Could try min..max. account which are the sweetspots. Weight is
                // the nr of positions this resolves.
                //
                // past readlen - kmerlen, or we could overcome
                if !self.n_stretch.is_set() && self.scp.is_past_contig() {
                    self.finalize_coding();
                }
                self.n_stretch.add_assign(1_u64);
                return None;
            }
        };
        if self.n_stretch.is_set() {
            self.finalize_n_stretch();
        }
        self.goffs.incr();
        //TODO: make sequence storage optional, to a separate file.
        self.ks.store_b2(&b2);
        if let Some(optimum_xmer) = self.scp.updated_median_xmer(&b2) {
            // FIXME: put the other strand in the idx top bits (saves reverse complementing).
            // Also, we need more alternative XmerLoc types for distinction along with traits !!
            // EDIT: not sure if the xor with revcmp really adds entropy.
            return self.is_recurrent(optimum_xmer);
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
                dbg_print!("sending {} to thread {}", mark, thread_index);
                if let Err(e) = self.tx[thread_index].send(mark) {
                    eprintln!("An sending to closed channel error here may mena that xmerhasher thread errored");
                    return Err(e.into());
                }
            }
        }
        dbg_print!("At end, processing last:");
        if self.n_stretch.is_set() {
            self.finalize_n_stretch();
        } else {
            self.finalize_coding();
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
}

impl<'a> XmerHash for KmerIter<'a> {
    fn get_kmerlen(&self) -> usize {
        usize::from(self.ks.get_bitlen()) / 2
    }
    fn get_overbit(&self) -> usize {
        let kmerlen = self.get_kmerlen();
        let dna_topb2_shift = (u32::try_from(kmerlen).unwrap() << 1) - 2;
        1_usize
            .checked_shl(dna_topb2_shift + 1)
            .expect("k-mer shift")
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::multi_thread;
    use crate::kmerconst::KmerConst;
    use crate::marker::fasta::Reader;
    use crate::new_types::extended_position::ExtPosEtc;
    use crate::new_types::position::BasePos;
    use anyhow::Result;
    use std::io;

    const GENOME_SIZE: u64 = 250;
    const READLEN: u16 = 6;

    fn read_from_string(s: &str) -> &[u8] {
        s.as_bytes()
    }

    fn process<'a>(ks: &'a mut KmerStore, kc: KmerConst, seq: &[u8]) -> Result<()> {
        let record = [b">test\n", seq].concat();
        let buf_reader = io::BufReader::new(&record[..]);
        let reader = Reader::new(buf_reader);
        multi_thread(ks, kc, reader, 2)
    }

    #[test]
    fn test_16n() -> Result<()> {
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, kc, b"NNNNNNNNNNNNNNNN")?;

        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, Position::default());
        dbg_assert_eq!(ks.contig[0].genomic, 16.into());
        Ok(())
    }
    #[test]
    fn test_1n() -> Result<()> {
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, kc, b"N")?;

        dbg_assert_eq!(ks.contig.len(), 1);
        dbg_assert_eq!(ks.contig[0].twobit, Position::default());
        dbg_assert_eq!(ks.contig[0].genomic, 1.into());
        Ok(())
    }
    #[test]
    fn test_1n1c1n() -> Result<()> {
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, kc, b"NCN")?;

        dbg_assert_eq!(ks.contig.len(), 2);
        // the position of the contig in ks.b2:
        dbg_assert_eq!(ks.contig[0].twobit, Position::default());
        // the position of the contig in genomic coordinates:
        dbg_assert_eq!(ks.contig[0].genomic, 1.into());

        // The last contig is a dummy, to indicate the end, again b2 position and genomic.
        dbg_assert_eq!(ks.contig[1].twobit, Position::from_basepos(1_u64));
        dbg_assert_eq!(ks.contig[1].genomic, 3.into());
        Ok(())
    }
    #[test]
    fn test_17c() -> Result<()> {
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let k = kc.kmerlen as u64;
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        process(&mut ks, kc, b"CCCCCCCCCCCCCCCCC")?;
        dbg_assert_eq!(ks.kmp.len(), 128);
        let mut first_pos = ExtPosEtc::from_basepos(k).get_rc();
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
    fn test_1n18c1n() -> Result<()> {
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let k = kc.kmerlen as u64;
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        {
            process(&mut ks, kc, b"NCCCCCCCCCCCCCCCCCCN")?;
        }
        let mut first_pos = ExtPosEtc::from_basepos(k).get_rc();
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
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let k = kc.kmerlen as u64;
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        {
            process(&mut ks, kc, b"NCCCCCCCCCCCCCCCC")?;
        }
        let mut seen = 0;
        let mut first_pos = ExtPosEtc::from_basepos(k).get_rc();
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
        let kc = KmerConst::new(GENOME_SIZE, READLEN, 0);
        let mut ks = KmerStore::new(kc.bitlen, 10_000, 0)?;
        {
            process(&mut ks, kc, b"ATATATATATATATATAT")?;
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
