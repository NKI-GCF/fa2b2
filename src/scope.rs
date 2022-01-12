use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{position::Position, twobit::TwoBit, xmer::Xmer};
use crate::rdbg::STAT_DB;
use anyhow::Result;

pub trait Scope {
    fn get_kc(&self) -> &KmerConst;
    fn get_d(&self, i: usize) -> &Xmer;
    fn update(&mut self) -> bool;
    fn pick_mark(&mut self) -> usize;
    fn dist_if_repetitive(
        &self,
        ks: &KmerStore,
        stored_p: ExtPosEtc,
        min_p: ExtPosEtc,
    ) -> Option<Position>;
    fn set_mark(&mut self, idx: usize, p: ExtPosEtc);
    fn increment(&mut self, b2: TwoBit);
}

pub trait WritingScope: Scope {
    fn is_repetitive(&self) -> bool;
    fn set_period(&mut self, period: Position);
    fn unset_period(&mut self);
    fn try_store_mark(
        &mut self,
        ks: &mut KmerStore,
        min_idx: usize,
        min_p: ExtPosEtc,
    ) -> Result<bool> {
        if self.is_repetitive() && ks.kmp[min_idx].is_set() {
            ks.kmp[min_idx].set_repetitive();
            return Ok(true);
        }
        let old_stored_p = ks.kmp[min_idx];

        if ks.kmp[min_idx].is_zero() {
            ks.set_kmp(min_idx, min_p);
            return Ok(true);
        }
        if old_stored_p.is_replaceable_by(min_p) {
            if old_stored_p.is_set_and_not(min_p) {
                // collision with a different baseindex, it had lesser extension.
                ks.set_kmp(min_idx, min_p);
                dbg_print!("[{:x}] -> {:?} (?)", min_idx, old_stored_p);
                while let Some((store_idx, store_p)) =
                    self.get_kc().get_next_xmer(min_idx, old_stored_p)
                {
                    if self.try_store_mark(ks, store_idx, store_p)? {
                        break;
                    }
                }
            }
            // .. else set and already min_p. Then leave bit states.
            return Ok(true);
        }
        if old_stored_p.extension() == min_p.extension() {
            // If a kmer occurs multiple times within an extending readlength (repetition),
            // only the first gets a position. During mapping this should be kept in mind.
            if let Some(dist) = self.dist_if_repetitive(ks, old_stored_p, min_p) {
                self.set_period(dist);
                // TODO: repetitive should be moved to higher extensions.
                ks.kmp[min_idx].set_repetitive();
                return Ok(true);
            }
            ks.kmp[min_idx].mark_more_recurs_upseq();
        }
        // collision with a different baseindex, it had greater extension.
        // the current baseindex will be extended and tried again.
        dbg_print!("not replacable, extend..");
        Ok(false)
    }

    fn store_mark(&mut self, ks: &mut KmerStore, i: usize) -> Result<()> {
        let mut min_idx = 0;
        let mut min_p = ExtPosEtc::zero();
        for x in 0..self.get_kc().get_ext_max() {
            (min_idx, min_p) = self.get_d(i).get_hash_and_p(x);
            if self.try_store_mark(ks, min_idx, min_p)? {
                break;
            }
        }
        self.set_mark(min_idx, min_p);
        Ok(())
    }
}
