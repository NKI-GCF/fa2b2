use crate::new_types::{
    extension::Extension,
    position::{BasePos, Position},
};
use crate::rdbg::STAT_DB;
use derive_more::Sub;
use serde::{Deserialize, Serialize};
use std::clone::Clone;
use std::ops::Add;
use std::{cmp, fmt};

const POS_SHIFT: u32 = 4;
const EXT_SHIFT: u32 = 56;
const ORI_MASK: u64 = 0x0000_0000_0000_0001;
const REP_MASK: u64 = 0x0000_0000_0000_0002;
const DUP_MASK: u64 = 0x0000_0000_0000_0004;
// TODO: indicate this position has annotation
const _INFO_MASK: u64 = 0x0000_0000_0000_0008;
const POS_MASK: u64 = 0x00FF_FFFF_FFFF_FFF0;
const EXT_MASK: u64 = 0xFF00_0000_0000_0000;

// FIXME u64 moet hier ExtPosEtc worden !!
impl Add<Position> for u64 {
    type Output = u64;
    fn add(self, other: Position) -> u64 {
        self + other.as_u64()
    }
}

#[derive(Sub, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
pub struct ExtPosEtc(u64);

impl ExtPosEtc {
    pub fn set(&mut self, p: ExtPosEtc) {
        *self = p;
    }
    pub fn as_u64(&self) -> u64 {
        self.0
    }
    pub fn basepos_to_pos(&self) -> Position {
        //FIXME: deprecate. as_pos() and pos() are confusing.
        Position::from(BasePos::from(*self))
    }
    pub fn unshift_pos(&self) -> u64 {
        //TODO: BasePos
        BasePos::from(self.pos()).as_u64()
    }
    pub fn pos(&self) -> Position {
        // until we have from::ExtPosEtc for
        Position::from(self.as_basepos())
    }
    pub fn as_basepos(&self) -> BasePos {
        //TODO: replace with BasePos::from()
        BasePos::from(*self)
    }
    pub fn zero() -> Self {
        ExtPosEtc(0x0)
    }
    pub fn clear(&mut self) {
        *self = ExtPosEtc::zero();
    }
    pub fn is_set(&self) -> bool {
        self.pos() != Position::zero()
    }
    pub fn is_zero(&self) -> bool {
        self.0 == 0
    }
    pub fn set_ori(&mut self, ori: bool) {
        self.0 ^= (self.0 ^ if ori { 1 } else { 0 }) & ORI_MASK
    }
    pub fn get_ori(&self) -> bool {
        self.0 & ORI_MASK != 0
    }
    pub fn incr_pos(&mut self) {
        self.0 += 1_u64.checked_shl(POS_SHIFT).expect("incr_pos shft");
    }
    pub fn decr_pos(&mut self) {
        dbg_assert!(self.is_set());
        self.0 -= 1_u64.checked_shl(POS_SHIFT).expect("decr_pos shft");
    }
    pub fn extension(&self) -> Extension {
        Extension::from(*self)
    }
    pub fn x(&self) -> usize {
        usize::from(self.extension())
    }
    pub fn same_ori(&self, p: ExtPosEtc) -> bool {
        self.get_ori() == p.get_ori()
    }
    pub fn blacklist(&mut self) {
        if self.is_set() {
            self.0 &= EXT_MASK;
            self.extend();
        }
    }
    pub fn extend(&mut self) {
        dbg_assert!(self.is_set());
        self.0 += 1_u64.checked_shl(EXT_SHIFT).expect("extend shft");
    }
    // preserves orientation, but not replication and duplication bits
    pub fn set_extension(&mut self, x: usize) {
        dbg_assert!(self.is_set());
        self.0 &= POS_MASK | ORI_MASK;
        self.0 |= Extension::from(x).as_u64();
    }
    pub fn clear_extension(&mut self) {
        self.0 &= !EXT_MASK;
    }
    pub fn is_same(&self, other: ExtPosEtc) -> bool {
        *self == other
    }
    /// Note: does not include ori, dup and rep bits.
    pub fn pos_with_ext(&self, x: usize) -> ExtPosEtc {
        ExtPosEtc(Extension::from(x).as_u64() | self.pos().as_u64())
    }
    pub fn set_dup(&mut self) {
        dbg_assert!(self.is_set());
        self.0 |= DUP_MASK;
    }
    pub fn is_dup(&self) -> bool {
        self.0 & DUP_MASK != 0
    }
    pub fn set_repetitive(&mut self) {
        dbg_assert!(self.is_set());
        self.0 |= REP_MASK;
    }
    pub fn is_repetitive(&self) -> bool {
        self.0 & REP_MASK != 0
    }
    pub fn rep_dup_masked(&self) -> ExtPosEtc {
        ExtPosEtc(self.0 & !(DUP_MASK | REP_MASK))
    }
    pub fn is_replaceable_by(&self, new_entry: ExtPosEtc) -> bool {
        // only extension bits means blacklisting, except for extension 0. pos is always > kmerlen
        new_entry.extension().as_u64() > self.rep_dup_masked().as_u64() // TODO: count down extension?
            || (new_entry.extension() == self.extension() && new_entry.pos() <= self.pos())
    }
    pub fn is_set_and_not(&self, other: ExtPosEtc) -> bool {
        self.is_set() && !self.is_same(other)
    }
    pub fn same_pos_and_ext(&self, new_entry: ExtPosEtc) -> bool {
        (self.0 ^ new_entry.0) & (EXT_MASK | POS_MASK) == 0
    }
    pub fn has_samepos(&self, other: ExtPosEtc) -> bool {
        self.pos() == other.pos() && {
            // XXX: may want to remove this later if orientation doesn't matter
            dbg_assert!(self.same_ori(other), "{:?}, {:?}", self, other);
            true
        }
    }
    /* does not work because extension is currently limited by readlen.
     * could maybe instead also shift base kmer or bits
     * fn ext_max() -> Self {
        0xFF
    }*/
}

impl From<Extension> for ExtPosEtc {
    fn from(e: Extension) -> ExtPosEtc {
        ExtPosEtc(e.as_u64())
    }
}

impl From<Position> for ExtPosEtc {
    fn from(pos: Position) -> ExtPosEtc {
        ExtPosEtc(pos.as_u64())
    }
}

impl From<BasePos> for ExtPosEtc {
    fn from(b: BasePos) -> ExtPosEtc {
        ExtPosEtc::from(Position::from(b))
    }
}

impl From<(Extension, Position)> for ExtPosEtc {
    fn from(ep: (Extension, Position)) -> ExtPosEtc {
        ExtPosEtc(ep.0.as_u64() | ep.1.as_u64())
    }
}

impl fmt::Debug for ExtPosEtc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:#x}", self.0)
    }
}

#[derive(new, Clone, PartialEq, Eq)]
pub struct KmerLoc {
    idx: usize,
    pub p: ExtPosEtc,
}
impl KmerLoc {
    pub fn get(&self) -> Option<(usize, ExtPosEtc)> {
        if self.is_set() {
            Some((self.idx, self.p))
        } else {
            None
        }
    }
    pub fn get_idx(&self) -> usize {
        self.idx
    }
    pub fn reset(&mut self) {
        self.idx = usize::max_value();
        self.p.clear();
    }

    pub fn is_set(&self) -> bool {
        self.idx != usize::max_value()
    }

    pub fn set(&mut self, idx: usize, p: ExtPosEtc, x: usize) {
        // during rebuilding the strange case occurs that mark is not set, but p is (extension)
        let self_p_extension = self.p.extension();
        let p_extension = p.extension();
        dbg_assert!(self.is_set() || self.p.is_zero() || self_p_extension == p_extension);
        self.idx = idx;
        self.p = p.rep_dup_masked();
        self.p.set_extension(x);
    }
}

impl PartialOrd for KmerLoc {
    fn partial_cmp(&self, other: &KmerLoc) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for KmerLoc {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        match self.idx.cmp(&other.idx) {
            /*FIXME: make extenion count down, etc to simplify to this:
            cmp::Ordering::Equal => match self.p.cmp(&other.p) {
                cmp::Ordering::Equal => panic!(),
                x => x,
            },*/
            cmp::Ordering::Equal => match other.p.extension().cmp(&self.p.extension()) {
                cmp::Ordering::Equal => match self.p.pos().cmp(&other.p.pos()) {
                    cmp::Ordering::Equal => panic!(),
                    x => x,
                },
                x => x,
            },
            x => x,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    impl KmerLoc {
        fn next(&mut self, ori: u64, is_template: bool) {
            if is_template {
                self.p.incr_pos()
            } else {
                self.p.decr_pos()
            }
            self.p.set_ori(ori & 1 != 0);
        }
    }
    #[test]
    fn forward() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, ExtPosEtc::from(BasePos::from(50_u64)));
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, true);
            dbg_assert_eq!(ori, if kl.p.get_ori() { 1 } else { 0 });
        }
        dbg_assert_eq!(
            kl.p.as_u64() & !1,
            ExtPosEtc::from(BasePos::from(50_u64 + pick)).as_u64()
        );
    }
    #[test]
    fn reverse() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, ExtPosEtc::from(BasePos::from(50_u64)));
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, false);
            dbg_assert_eq!(ori, if kl.p.get_ori() { 1 } else { 0 });
        }
        dbg_assert_eq!(
            kl.p.as_u64() & !1,
            ExtPosEtc::from(BasePos::from(50_u64 - pick)).as_u64()
        );
    }
}
