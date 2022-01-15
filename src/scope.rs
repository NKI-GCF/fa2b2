use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{
    position::Position,
    twobit::{ThreeBit, TwoBit},
    xmer::Xmer,
};
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use anyhow::Result;

pub trait Scope {
    fn get_pos(&self) -> Position;
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
    fn set_mark(&mut self, mark: &XmerLoc);
    fn increment(&mut self, b2: TwoBit);
    fn ascii_to_b3(&self, b: &u8) -> ThreeBit {
        ThreeBit::from((self.get_pos(), *b))
    }
}

pub trait WritingScope: Scope {
    fn is_repetitive(&self) -> bool;
    fn set_period(&mut self, period: Position);
    fn unset_period(&mut self);

    fn try_store_mark(&mut self, ks: &mut KmerStore, mark: &mut XmerLoc) -> Result<bool> {
        if self.is_repetitive() && ks.kmp[mark.idx].is_set() {
            ks.kmp[mark.idx].set_repetitive();
            return Ok(false);
        }
        let old_stored_p = ks.kmp[mark.idx];

        if old_stored_p.is_zero() {
            ks.set_kmp(&mark);
            return Ok(false);
        }
        if old_stored_p.is_replaceable_by(mark.p) {
            if old_stored_p.pos() != mark.p.pos() {
                ks.set_kmp(&mark);
                if old_stored_p.x() == mark.p.x() {
                    // same extension means same base k-mer origin. this is a duplicate.
                    ks.kmp[mark.idx].mark_more_recurs_upseq();
                }
                mark.p = old_stored_p;
                dbg_print!("[{:x}] -> {:?} (?)", mark.idx, mark.p);
                return Ok(true);
            }
            // .. else set and already mark.p. Then leave bit states.
            return Ok(false);
        }
        if old_stored_p.extension() == mark.p.extension() {
            // this must be the same k-mer origin, meaning identical k-mer sequence.

            // If a kmer occurs multiple times within an extending readlength (repetition),
            // only the first gets a position. During mapping this should be kept in mind.
            if let Some(dist) = self.dist_if_repetitive(ks, old_stored_p, mark.p) {
                self.set_period(dist);
                // TODO: repetitive should be moved to higher extensions.
                ks.kmp[mark.idx].set_repetitive();
                return Ok(false);
            }
            ks.kmp[mark.idx].mark_more_recurs_upseq();
        }
        // collision with a different baseindex, it had greater extension.
        // the current baseindex will be extended and tried again.
        dbg_print!("not replacable, extend..");
        Ok(true)
    }

    fn store_mark(&mut self, ks: &mut KmerStore, i: usize) -> Result<()> {
        let mut mark = self.get_d(i).get_hash_and_p(self.get_kc(), 0);
        self.set_mark(&mark);
        let orig_pos = mark.p.pos();
        // this seems to be hotlooping
        while self.try_store_mark(ks, &mut mark)? {
            if mark.p.pos() != orig_pos {
                if self.get_kc().extend_xmer(&mut mark).is_ok() {
                    // extending some pase baseidx. TODO: if frequently the same recurs,
                    // it might be worthwhile to store the reverse complement in a temp
                    // we should not store a past index in self.mark.p !!
                } else {
                    dbg_print!("couldn't extend: [{}] {:?} ..?", mark.idx, mark.p);
                    break;
                }
            } else {
                if mark.p.extend().is_err() {
                    break;
                }
                mark = self.get_d(i).get_hash_and_p(self.get_kc(), mark.p.x());
                self.set_mark(&mark);
            }
        }
        Ok(())
    }
}
