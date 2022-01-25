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

// lowering this below 8 may cause some test failures
//const EXT_BITS: usize = 8;

// FIXME: ext_bits hangt af van de kmersize.
// FIXME 2: zonder ori bit is de extension een bit langer, dat maakt vanalles eenvoudiger.
pub(crate) const EXT_BITS: usize = 16;

pub(crate) const EXT_MAX: usize = 1 << EXT_BITS;

pub(super) const EXT_SHIFT: u32 = (size_of::<u64>() * 8 - EXT_BITS) as u32;
pub(crate) const EXT_MASK: u64 = !0x0 ^ ((1 << EXT_SHIFT) - 1);
const REP_SHIFT: u32 = EXT_SHIFT - 1;
const DUP_SHIFT: u32 = EXT_SHIFT - 2;
const ETC_SHIFT: u32 = EXT_SHIFT - 3;
pub(crate) const REPETITIVE: u64 = 1 << REP_SHIFT;
pub(crate) const DUPLICATE: u64 = 1 << DUP_SHIFT;
pub(crate) const ETCETERA: u64 = 1 << ETC_SHIFT;
pub(super) const POS_SHIFT: u32 = 3;
pub(super) const POS_MASK: u64 = (ETCETERA - 1) ^ ((1 << POS_SHIFT) - 1);

const ORI_MASK: u64 = 0x1;
// TODO: indicate this position has annotation

#[derive(
    Display,
    DebugCustom,
    Sub,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
    Hash,
    Default,
)]
#[display(fmt = "{:#x}", _0)]
#[debug(fmt = "{:#x}", _0)]
// needs to be pub or "can't leak crate-private type"
pub struct ExtPosEtc(u64);

#[derive(Default, Clone)]
pub(crate) struct MiniExtPosOri(u32);
impl MiniExtPosOri {}

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
    /// a position is shifted, but has no other thant the POS_MASK bits set
    pub(crate) fn pos(&self) -> Position {
        Position::from(*self)
    }
    pub(crate) fn from_basepos<T>(value: T) -> ExtPosEtc
    where
        T: Into<u64>,
    {
        ExtPosEtc::from((Extension::default(), Position::from_basepos(value)))
    }
    #[deprecated]
    pub(crate) fn basepos(&self) -> BasePos {
        BasePos::from(self.pos())
    }
    pub(crate) fn as_basepos(&self) -> BasePos {
        BasePos::from(*self)
    }
    pub(crate) fn clear(&mut self) {
        *self = ExtPosEtc::default();
    }
    pub(crate) fn is_set(&self) -> bool {
        self.pos() != Position::default()
    }
    #[inline(always)]
    pub(crate) fn is_zero(&self) -> bool {
        self.0 == 0
    }
    pub(crate) fn set_ori(&mut self, ori: bool) {
        self.0 ^= (self.0 ^ if ori { 1 } else { 0 }) & ORI_MASK
    }
    #[deprecated] // use !is_template
    pub(crate) fn get_ori(&self) -> bool {
        self.0 & ORI_MASK != 0
    }
    pub(crate) fn is_template(&self) -> bool {
        self.0 & ORI_MASK == 0
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
        self.is_template() == p.is_template()
    }
    pub(crate) fn get_rc(&self) -> ExtPosEtc {
        ExtPosEtc(self.0 ^ ORI_MASK)
    }
    pub(crate) fn blacklist(&mut self) {
        if self.is_set() {
            self.0 &= EXT_MASK;
            self.extend().unwrap();
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn no_mask_overlaps() {
        dbg_assert_eq!(EXT_MASK & POS_MASK, 0);
        dbg_assert_eq!(EXT_MASK & ORI_MASK, 0);
        dbg_assert_eq!(EXT_MASK & REPETITIVE, 0);
        dbg_assert_eq!(EXT_MASK & DUPLICATE, 0);
        dbg_assert_eq!(EXT_MASK & ETCETERA, 0);

        dbg_assert_eq!(POS_MASK & ORI_MASK, 0);
        dbg_assert_eq!(POS_MASK & REPETITIVE, 0);
        dbg_assert_eq!(POS_MASK & DUPLICATE, 0);
        dbg_assert_eq!(POS_MASK & ETCETERA, 0);

        dbg_assert_eq!(ORI_MASK & REPETITIVE, 0);
        dbg_assert_eq!(ORI_MASK & DUPLICATE, 0);
        dbg_assert_eq!(ORI_MASK & ETCETERA, 0);

        dbg_assert_eq!(REPETITIVE & DUPLICATE, 0);
        dbg_assert_eq!(REPETITIVE & ETCETERA, 0);

        dbg_assert_eq!(DUPLICATE & ETCETERA, 0);
    }
}
