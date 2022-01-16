use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{
    position::Position,
    twobit::{ThreeBit, TwoBit},
};
use crate::xmer_location::XmerLoc;

pub trait Scope {
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
}
