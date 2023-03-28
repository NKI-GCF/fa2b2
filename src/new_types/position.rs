// (c) Roel Kluin, 2023, GPL v3

use crate::new_types::extended_position::{ExtPosEtc, POS_MASK, POS_SHIFT};
use crate::num::ToPrimitive;
use crate::rdbg::STAT_DB;
use anyhow::Error as AnyhowError;
use anyhow::{anyhow, ensure, Result};
use derive_more::{Add, AddAssign, Rem, Sub};
use serde::{Deserialize, Serialize};
use std::clone::Clone;
use std::convert::TryFrom;
use std::fmt;

// twobit shifts are bit positions 0, 2, 4 and 6 in a byte.
const TWOBIT_SHIFT: u32 = POS_SHIFT - 1;

const BASEPOS_MASK: u64 = POS_MASK >> POS_SHIFT;
// Two types are defined here.
//
/// a Position is an u64 with only bits in POS_MASK set. Upper and lower bits are 0.
/// The position is a contig-cumulative genmic offset, currected for any stretches of Ns.
#[derive(
    Add,
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
    Rem,
    Display,
    DebugCustom,
    Default,
)]
#[display(fmt = "{:#x}", _0)]
#[debug(fmt = "{:#x}", _0)]
pub struct Position(u64);

// only bits set for postion, but not shifted yet
/// a BasePos is a [`Position`] shifted right with POS_SHIFT. High bits are unset.
#[derive(
    Add, Sub, Clone, Copy, Serialize, Deserialize, AddAssign, Debug, Display, Default, PartialEq, Eq,
)]
pub struct BasePos(u64);

#[derive(Clone, Copy, Default)]
pub(crate) struct PosRange((Position, Position));

impl PosRange {
    pub(crate) fn lower(&self) -> Position {
        self.0 .0
    }
    pub(crate) fn upper(&self) -> Position {
        self.0 .1
    }
    pub(crate) fn has_in_range(&self, pos: Position) -> bool {
        pos >= self.0 .0 && pos < self.0 .1
    }
    pub(crate) fn bound_lower(&mut self, value: Position) {
        let inner = &mut self.0;
        if value > inner.0 {
            dbg_assert!(value < inner.1);
            inner.0 = value;
        }
    }
    pub(crate) fn bound_upper(&mut self, value: Position) {
        let inner = &mut self.0;
        if value < inner.1 {
            dbg_assert!(value > inner.0);
            inner.1 = value;
        }
    }
}

impl TryFrom<(Position, Position)> for PosRange {
    type Error = AnyhowError;
    fn try_from(rng: (Position, Position)) -> Result<PosRange, Self::Error> {
        ensure!(rng.0 < rng.1);
        Ok(PosRange(rng))
    }
}

impl fmt::Display for PosRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}-{}",
            BasePos::from(self.0 .0),
            BasePos::from(self.0 .1)
        )
    }
}

impl Position {
    pub(crate) fn as_u64(&self) -> u64 {
        self.0
    }
    pub(crate) fn is_set(&self) -> bool {
        self.0 != 0
    }
    // before applying TWOBIT_SHIFT (4 places) I had this:
    // XXX: thread 'main' panicked at 'u32??: TryFromIntError(())', src/new_types/position.rs:49:39
    // TODO: check that distances are sensible.
    pub(crate) fn b2_shift(&self) -> u32 {
        let bit_pos = self.0.checked_shr(TWOBIT_SHIFT).unwrap() & 6;
        bit_pos.to_u32().unwrap()
    }

    // for a repetitive mark, return if on period for the distance between current pos (self) and stored
    pub(crate) fn get_if_mark_on_period(
        &self,
        stored_pos: Position,
        pd: Position,
    ) -> Option<BasePos> {
        let dist = if self.0 >= stored_pos.0 {
            self.0 - stored_pos.0
        } else {
            stored_pos.0 - self.0
        };
        if dist % pd.0 == 0 {
            Some(BasePos(dist >> POS_SHIFT))
        } else {
            None
        }
    }
    // increment the position one
    pub(crate) fn incr(&mut self) {
        self.0 += 1_u64.checked_shl(POS_SHIFT).expect("pos.incr shft");
        dbg_assert_eq!(self.0, self.0 & POS_MASK);
    }
    // increment the position one
    pub(crate) fn decr(&mut self) {
        dbg_assert_ne!(self.0, 0);
        self.0 -= 1_u64.checked_shl(POS_SHIFT).expect("pos.incr shft");
    }
    pub(crate) fn from_basepos<T>(value: T) -> Position
    where
        T: Into<u64>,
    {
        Position::from(BasePos::from(value.into()))
    }
}

impl From<BasePos> for Position {
    fn from(base_pos: BasePos) -> Position {
        Position(u64::from(base_pos) << POS_SHIFT)
    }
}

impl From<ExtPosEtc> for Position {
    fn from(p: ExtPosEtc) -> Position {
        Position(p.as_u64() & POS_MASK)
    }
}

/////////////// BasePos \\\\\\\\\\\\\\\\\\

impl BasePos {
    pub(crate) fn as_u64(&self) -> u64 {
        self.0
    }
    pub(crate) fn is_set(&self) -> bool {
        self.0 != 0
    }
    pub(crate) fn as_usize(&self) -> usize {
        usize::try_from(self.0).unwrap()
    }
    pub(crate) fn incr(&mut self) {
        self.0 += 1;
    }
    #[must_use = "this value should be used"]
    pub(crate) fn add<T>(&mut self, add: T) -> BasePos
    where
        T: Into<u64>,
    {
        let add = add.into();
        dbg_assert!(self.0 + add <= BASEPOS_MASK);
        BasePos(self.0 + add)
    }
    pub(crate) fn add_assign<T>(&mut self, add: T)
    where
        T: Into<u64>,
    {
        *self = self.add(add);
    }
}

// converts into base pos.
impl From<Position> for BasePos {
    fn from(pos: Position) -> BasePos {
        BasePos(pos.as_u64().checked_shr(POS_SHIFT).unwrap())
    }
}

// converts into base pos.
impl From<ExtPosEtc> for BasePos {
    fn from(val: ExtPosEtc) -> BasePos {
        BasePos::from(Position::from(val))
    }
}

impl From<u64> for BasePos {
    fn from(base_pos: u64) -> BasePos {
        dbg_assert!(base_pos <= BASEPOS_MASK, "bleeds beyond position!");
        BasePos(base_pos)
    }
}

macro_rules! implement_try_froms_for_basepos { ($($ty:ty),*) => ($(
    impl TryFrom<$ty> for BasePos
    where $ty: TryFrom<u64> {
        type Error = AnyhowError;
        fn try_from(base_pos: $ty) -> Result<BasePos> {
            u64::try_from(base_pos).map(BasePos::from).map_err(|e| anyhow!("{}", e))
        }
    }
    )*)
}
implement_try_froms_for_basepos!(usize, i64, u32);

impl From<BasePos> for u64 {
    fn from(base_pos: BasePos) -> u64 {
        base_pos.0
    }
}
