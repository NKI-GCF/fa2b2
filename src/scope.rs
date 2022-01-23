use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{position::Position, twobit::TwoBit};

pub trait Scope {
    fn update(&mut self);
    #[deprecated]
    fn dist_if_repetitive(
        &self,
        ks: &KmerStore,
        stored_p: ExtPosEtc,
        min_p: ExtPosEtc,
    ) -> Option<Position>;
    fn set_mark(&mut self, i: usize);
    fn increment(&mut self, b2: TwoBit);
}
