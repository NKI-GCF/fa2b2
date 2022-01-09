use crate::new_types::{
    extension::Extension,
    position::{BasePos, Position},
};
use crate::rdbg::STAT_DB;
use std::clone::Clone;
use std::cmp;
use std::ops::Add;
//use std::convert::TryFrom;
//use serde::{Deserialize, Serialize};

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

// FIXME, should become a newtype rather than a trait.
pub trait ExtPosEtc: Clone {
    fn set(&mut self, p: u64);
    fn get(&self) -> u64;
    fn pos(&self) -> Position;
    fn as_basepos(&self) -> BasePos;
    fn basepos_to_pos(&self) -> Position;
    fn unshift_pos(&self) -> u64;
    fn zero() -> Self;
    fn clear(&mut self);
    fn is_set(&self) -> bool;
    fn is_zero(&self) -> bool;
    fn set_ori(&mut self, p: u64);
    fn get_ori(&self) -> u64;
    fn incr_pos(&mut self);
    fn decr_pos(&mut self);
    fn extension(&self) -> Extension;
    fn x(&self) -> usize;
    fn same_ori(&self, p: u64) -> bool;
    fn blacklist(&mut self);
    fn extend(&mut self);
    fn set_extension(&mut self, x: u64);
    fn clear_extension(&mut self);
    fn is_same(&self, other: u64) -> bool;
    fn pos_with_ext(&self, x: usize) -> u64;
    fn set_dup(&mut self);
    fn is_dup(&self) -> bool;
    fn set_repetitive(&mut self);
    fn is_repetitive(&self) -> bool;
    fn rep_dup_masked(&self) -> u64;
    fn is_replaceable_by(&self, new_entry: u64) -> bool;
    fn is_set_and_not(&self, other: u64) -> bool;
    fn same_pos_and_ext(&self, new_entry: u64) -> bool;
    fn has_samepos(&self, other: u64) -> bool;
    //fn ext_max() -> Self;
}

impl ExtPosEtc for u64 {
    fn set(&mut self, p: u64) {
        *self = p;
    }
    fn get(&self) -> u64 {
        dbg_assert!(self.is_set());
        *self
    }
    fn basepos_to_pos(&self) -> Position {
        //FIXME: deprecate. as_pos() and pos() are confusing.
        Position::from(BasePos::from(*self))
    }
    fn unshift_pos(&self) -> u64 {
        //TODO: BasePos
        //FIXME: change so dup, rep etc. is in first nibble.
        BasePos::from(self.pos()).as_u64()
    }
    fn pos(&self) -> Position {
        // until we have from::ExtPosEtc for
        Position::from(self.as_basepos())
    }
    fn as_basepos(&self) -> BasePos {
        BasePos::from(self.checked_shr(POS_SHIFT).expect("shr pos shft"))
    }
    fn zero() -> Self {
        0x0
    }
    fn clear(&mut self) {
        *self = ExtPosEtc::zero();
    }
    fn is_set(&self) -> bool {
        self.pos() != Position::zero()
    }
    fn is_zero(&self) -> bool {
        self.pos() == Position::zero()
    }
    fn set_ori(&mut self, p: u64) {
        *self ^= (*self ^ p) & ORI_MASK
    }
    fn get_ori(&self) -> u64 {
        dbg_assert!(self.is_set());
        self & ORI_MASK
    }
    fn incr_pos(&mut self) {
        *self += 1_u64.checked_shl(POS_SHIFT).expect("incr_pos shft");
    }
    fn decr_pos(&mut self) {
        dbg_assert!(self.is_set());
        *self -= 1_u64.checked_shl(POS_SHIFT).expect("decr_pos shft");
    }
    // not all extensions may apply, it's dependent on genome size.
    fn extension(&self) -> Extension {
        // TODO: rewrite when ExtPosEtc is a new type
        //FIXME: change so that extension is counted down.

        Extension::from(((*self & EXT_MASK) >> EXT_SHIFT) as usize)
    }
    fn x(&self) -> usize {
        usize::from(self.extension())
    }
    fn same_ori(&self, p: u64) -> bool {
        dbg_assert!(self.is_set());
        self.get_ori() == (p & ORI_MASK)
    }
    fn blacklist(&mut self) {
        if self.is_set() {
            *self &= EXT_MASK;
            self.extend();
        }
    }
    fn extend(&mut self) {
        dbg_assert!(self.is_set());
        *self += 1_u64.checked_shl(EXT_SHIFT).expect("extend shft");
    }
    fn set_extension(&mut self, x: u64) {
        dbg_assert!(self.is_set());
        //XXX does not work with !EXT_MASK; (?)
        *self &= POS_MASK | ORI_MASK;
        *self |= x.checked_shl(EXT_SHIFT).expect("x beyond max extension");
    }
    fn clear_extension(&mut self) {
        dbg_assert!(self.is_set());
        *self &= !EXT_MASK;
    }
    fn is_same(&self, other: u64) -> bool {
        dbg_assert!(self.is_set());
        *self == other
    }
    /// Note: unsets dup and rep bits.
    fn pos_with_ext(&self, x: usize) -> u64 {
        self.pos().as_u64() | (x as u64).checked_shl(EXT_SHIFT).expect("x beyond max ext")
    }
    fn set_dup(&mut self) {
        dbg_assert!(self.is_set());
        *self |= DUP_MASK;
    }
    fn is_dup(&self) -> bool {
        dbg_assert!(self.is_set());
        *self & DUP_MASK != 0
    }
    fn set_repetitive(&mut self) {
        dbg_assert!(self.is_set());
        *self |= REP_MASK;
    }
    fn is_repetitive(&self) -> bool {
        dbg_assert!(self.is_set());
        *self & REP_MASK != 0
    }
    fn rep_dup_masked(&self) -> u64 {
        dbg_assert!(self.is_set());
        *self & !(DUP_MASK | REP_MASK)
    }
    fn is_replaceable_by(&self, new_entry: u64) -> bool {
        // only extension bits means blacklisting, except for extension 0. pos is always > kmerlen
        new_entry.extension().as_u64() > self.rep_dup_masked() // TODO: count down extension?
            || (new_entry.extension() == self.extension() && new_entry.pos() <= self.pos())
    }
    fn is_set_and_not(&self, other: u64) -> bool {
        self.is_set() && !self.is_same(other)
    }
    fn same_pos_and_ext(&self, new_entry: u64) -> bool {
        dbg_assert!(self.is_set());
        (*self ^ new_entry) & (EXT_MASK | POS_MASK) == 0
    }
    fn has_samepos(&self, other: u64) -> bool {
        self.pos() == other.pos() && {
            // XXX: may want to remove this later if orientation doesn't matter
            dbg_assert!(self.same_ori(other), "{:#x}, {:#x}", self, other);
            true
        }
    }
    /* does not work because extension is currently limited by readlen.
     * could maybe instead also shift base kmer or bits
     * fn ext_max() -> Self {
        0xFF
    }*/
}

#[derive(new, Clone, PartialEq, Eq)]
pub struct KmerLoc {
    idx: usize,
    pub p: u64,
}
impl KmerLoc {
    pub fn get(&self) -> (usize, u64) {
        (self.idx, self.p)
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

    pub fn set(&mut self, idx: usize, p: u64, x: usize) {
        // during rebuilding the strange case occurs that mark is not set, but p is (extension)
        let self_p_extension = self.p.extension();
        let p_extension = p.extension();
        dbg_assert!(self.is_set() || self.p.is_zero() || self_p_extension == p_extension);
        self.idx = idx;
        self.p = p.rep_dup_masked();
        self.p.set_extension(x as u64);
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
        fn next(&mut self, p: u64, is_template: bool) {
            if is_template {
                self.p.incr_pos()
            } else {
                self.p.decr_pos()
            }
            self.p.set_ori(p);
        }
    }
    #[test]
    fn forward() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 50.basepos_to_pos().as_u64());
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, true);
            dbg_assert_eq!(ori, kl.p & 1);
        }
        dbg_assert_eq!(kl.p & !1, (50 + pick).basepos_to_pos().as_u64());
    }
    #[test]
    fn reverse() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 50.basepos_to_pos().as_u64());
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, false);
            dbg_assert_eq!(ori, kl.p & 1);
        }
        dbg_assert_eq!(kl.p & !1, (50 - pick).basepos_to_pos().as_u64());
    }
}
