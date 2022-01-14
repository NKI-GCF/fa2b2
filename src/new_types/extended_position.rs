use crate::new_types::{
    extension::Extension,
    position::{BasePos, Position},
};
use crate::rdbg::STAT_DB;
use anyhow::{anyhow, Result};
use derive_more::Sub;
use serde::{Deserialize, Serialize};
use std::clone::Clone;
use std::mem::size_of;
use std::{cmp, fmt};

// lowering this below 8 may cause some test failures
const EXT_BITS: usize = 8;
//TODO: const EXT_BITS: u32 = 16;

pub const EXT_MAX: usize = 1 << EXT_BITS;

pub(super) const EXT_SHIFT: u32 = (size_of::<u64>() * 8 - EXT_BITS) as u32;
pub(crate) const EXT_MASK: u64 = !0x0 ^ ((1 << EXT_SHIFT) - 1);
const REP_SHIFT: u32 = EXT_SHIFT - 1;
pub(crate) const REPETITIVE: u64 = 1 << REP_SHIFT;
pub(super) const POS_SHIFT: u32 = 3;
pub(super) const POS_MASK: u64 = (REPETITIVE - 1) ^ ((1 << POS_SHIFT) - 1);

const ORI_MASK: u64 = 0x1;
pub(crate) const DUPLICATE: u64 = 0x2;
// TODO: indicate this position has annotation
const _INFO_MASK: u64 = 0x4;

#[derive(Sub, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
pub struct ExtPosEtc(u64);

/// The stored information per xmer: in u64 position, extension and flags:
/// Orientation, Duplicate, Replication(, HasInfo: TODO)
impl ExtPosEtc {
    pub(crate) fn set(&mut self, p: ExtPosEtc) {
        *self = p;
    }
    pub(crate) fn as_u64(&self) -> u64 {
        self.0
    }
    pub(crate) fn unshift_pos(&self) -> u64 {
        BasePos::from(self.pos()).as_u64()
    }
    pub(crate) fn pos(&self) -> Position {
        // until we have from::ExtPosEtc for
        Position::from(self.as_basepos())
    }
    pub(crate) fn as_basepos(&self) -> BasePos {
        BasePos::from(*self)
    }
    pub(crate) fn zero() -> Self {
        ExtPosEtc(0x0)
    }
    pub(crate) fn clear(&mut self) {
        *self = ExtPosEtc::zero();
    }
    pub(crate) fn is_set(&self) -> bool {
        self.pos() != Position::zero()
    }
    #[inline(always)]
    pub(crate) fn is_zero(&self) -> bool {
        self.0 == 0
    }
    pub(crate) fn set_ori(&mut self, ori: bool) {
        self.0 ^= (self.0 ^ if ori { 1 } else { 0 }) & ORI_MASK
    }
    pub(crate) fn get_ori(&self) -> bool {
        self.0 & ORI_MASK != 0
    }
    pub(crate) fn incr_pos(&mut self) {
        self.0 += 1_u64.checked_shl(POS_SHIFT).expect("incr_pos shft");
    }
    pub(crate) fn decr_pos(&mut self) {
        dbg_assert!(self.is_set());
        self.0 -= 1_u64.checked_shl(POS_SHIFT).expect("decr_pos shft");
    }
    pub(crate) fn extension(&self) -> Extension {
        Extension::from(*self)
    }
    pub(crate) fn x(&self) -> usize {
        usize::from(self.extension())
    }
    pub(crate) fn same_ori(&self, p: ExtPosEtc) -> bool {
        self.get_ori() == p.get_ori()
    }
    pub(crate) fn blacklist(&mut self) {
        if self.is_set() {
            self.0 &= EXT_MASK;
            self.extend();
        }
    }
    pub(crate) fn extend(&mut self) -> Result<()> {
        dbg_assert!(self.is_set());
        self.0 = self
            .0
            .checked_add(1_u64.checked_shl(EXT_SHIFT).expect("extend shft"))
            .ok_or_else(|| anyhow!("at max extension"))?;
        Ok(())
    }
    pub(crate) fn clear_extension(&mut self) {
        self.0 &= !EXT_MASK;
    }
    pub(crate) fn mark_more_recurs_upseq(&mut self) {
        dbg_assert!(self.is_set());
        self.0 |= DUPLICATE;
    }
    pub(crate) fn is_dup(&self) -> bool {
        self.0 & DUPLICATE != 0
    }
    pub(crate) fn is_last_on_ref(&self) -> bool {
        !self.is_dup()
    }
    pub(crate) fn set_repetitive(&mut self) {
        dbg_assert!(self.is_set());
        self.0 |= REPETITIVE;
    }
    pub(crate) fn is_repetitive(&self) -> bool {
        self.0 & REPETITIVE != 0
    }
    #[must_use]
    pub(crate) fn rep_dup_masked(&self) -> ExtPosEtc {
        ExtPosEtc(self.0 & !(DUPLICATE | REPETITIVE))
    }
    pub(crate) fn is_replaceable_by(&self, new_entry: ExtPosEtc) -> bool {
        // only extension bits means blacklisting, except for extension 0. pos is always > kmerlen
        new_entry.extension().as_u64() > self.rep_dup_masked().as_u64() // TODO: count down extension?
            || (new_entry.extension() == self.extension() && new_entry.pos() <= self.pos())
    }
    pub(crate) fn is_set_and_not(&self, other: ExtPosEtc) -> bool {
        self.is_set() && self.pos() != other.pos()
    }
    pub(crate) fn same_pos_and_ext(&self, new_entry: ExtPosEtc) -> bool {
        (self.0 ^ new_entry.0) & (EXT_MASK | POS_MASK) == 0
    }
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
    pub(crate) fn get(&self) -> Option<(usize, ExtPosEtc)> {
        if self.is_set() {
            Some((self.idx, self.p))
        } else {
            None
        }
    }
    pub(crate) fn get_idx(&self) -> usize {
        self.idx
    }
    pub(crate) fn reset(&mut self) {
        self.idx = usize::max_value();
        self.p.clear();
    }

    pub(crate) fn is_set(&self) -> bool {
        self.idx != usize::max_value()
    }

    pub(crate) fn set(&mut self, idx: usize, p: ExtPosEtc) {
        // during rebuilding the strange case occurs that mark is not set, but p is (extension)
        let self_p_extension = self.p.extension();
        let p_extension = p.extension();
        dbg_assert!(self.is_set() || self.p.is_zero() || self_p_extension == p_extension);
        self.idx = idx;
        self.p = p;
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
    fn no_mask_overlaps() {
        dbg_assert_eq!(EXT_MASK & POS_MASK, 0);
        dbg_assert_eq!(EXT_MASK & ORI_MASK, 0);
        dbg_assert_eq!(EXT_MASK & REPETITIVE, 0);
        dbg_assert_eq!(EXT_MASK & DUPLICATE, 0);
        dbg_assert_eq!(EXT_MASK & _INFO_MASK, 0);

        dbg_assert_eq!(POS_MASK & ORI_MASK, 0);
        dbg_assert_eq!(POS_MASK & REPETITIVE, 0);
        dbg_assert_eq!(POS_MASK & DUPLICATE, 0);
        dbg_assert_eq!(POS_MASK & _INFO_MASK, 0);

        dbg_assert_eq!(ORI_MASK & REPETITIVE, 0);
        dbg_assert_eq!(ORI_MASK & DUPLICATE, 0);
        dbg_assert_eq!(ORI_MASK & _INFO_MASK, 0);

        dbg_assert_eq!(REPETITIVE & DUPLICATE, 0);
        dbg_assert_eq!(REPETITIVE & _INFO_MASK, 0);

        dbg_assert_eq!(DUPLICATE & _INFO_MASK, 0);
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
