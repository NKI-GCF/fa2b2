use crate::kmerconst::KmerConst;
use crate::kmerconst::XmerHash;
use crate::new_types::extended_position::ExtPosEtc;
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use anyhow::Result;

pub struct XmerHashExtensionThread {
    kmp: Vec<ExtPosEtc>, // position + strand per xmer handled in this thread.
    thread_nr: usize,
    thread_max: usize,
    kmerlen: usize,
    overbit: usize,
}

impl XmerHashExtensionThread {
    pub(crate) fn new(kc: &KmerConst, thread_nr: u32, thread_max: u32) -> Self {
        dbg_assert!(thread_max.is_power_of_two());
        let shl = kc.bitlen - 2 - usize::try_from(thread_max.trailing_zeros()).unwrap();
        XmerHashExtensionThread {
            kmp: vec![ExtPosEtc::default(); 1 << shl],
            thread_nr: usize::try_from(thread_nr).unwrap(),
            thread_max: usize::try_from(thread_max).unwrap(),
            kmerlen: kc.kmerlen,
            overbit: kc.overbit,
        }
    }
    fn store_mark_and_extend(&mut self, mark: &mut XmerLoc) -> Result<()> {
        while self.try_store_mark(mark)? {
            let orig_pos = mark.p.pos();
            if mark.p.pos() != orig_pos {
                if self.extend_xmer(mark).is_ok() {
                    // extending some pase baseidx. TODO: if frequently the same recurs,
                    // it might be worthwhile to store the reverse complement in a temp
                    // we should not store a past index in self.mark.p !!
                } else {
                    dbg_print!("couldn't extend: {} ..?", mark);
                    break;
                }
            } else {
                if self.extend_xmer(mark).is_err() {
                    break;
                }
                if !self.is_handled_in_thread(&mark) {
                    // XXX send to another thread
                    break;
                }
            }
        }
        Ok(())
    }
    /// the top bits of the xmerhash
    pub(crate) fn is_handled_in_thread(&self, mark: &XmerLoc) -> bool {
        let kmer = self.unhash_and_uncompress_to_kmer(mark.idx, mark.p.x());
        kmer & (self.thread_max - 1) == self.thread_nr
    }
    fn set_kmp(&mut self, mark: &XmerLoc) {
        self.kmp[mark.idx].set(mark.p);
    }
    fn try_store_mark(&mut self, mark: &mut XmerLoc) -> Result<bool> {
        let old_stored_p = self.kmp[mark.idx];

        if old_stored_p.is_zero() {
            self.set_kmp(&mark);
            return Ok(false);
        }
        if old_stored_p.is_replaceable_by(mark.p) {
            if old_stored_p.pos() == mark.p.pos() {
                // set and already mark.p. Leave the bit states.
                return Ok(false);
            }
            self.set_kmp(&mark);
            if old_stored_p.x() == mark.p.x() {
                // same extension means same base k-mer origin. this is a duplicate.
                self.kmp[mark.idx].mark_more_recurs_upseq();
            }
            dbg_print!("{} -> ?", mark);
            mark.p = old_stored_p; // to be extended next
        } else if old_stored_p.extension() == mark.p.extension() {
            // Note: same extension and hash means same k-mer origin: identical k-mer sequence.
            self.kmp[mark.idx].mark_more_recurs_upseq();
        }
        // collision between hashes, the one in mark.p will be extended and tried again.
        Ok(true)
    }
}

impl XmerHash for XmerHashExtensionThread {
    fn get_kmerlen(&self) -> usize {
        self.kmerlen
    }
    fn get_overbit(&self) -> usize {
        self.overbit
    }
}
