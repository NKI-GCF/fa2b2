use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{
    position::Position,
    twobit::{PosB3, ThreeBit, TwoBit},
};
use crate::xmer_location::XmerLoc;

pub trait Scope {
    fn get_pos(&self) -> Position;
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
    fn ascii_to_b3(&self, b: &u8) -> PosB3 {
        ThreeBit::get_pos_b3(self.get_pos(), *b)
    }
}
