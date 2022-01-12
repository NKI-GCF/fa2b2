use crate::new_types::extended_position::{ExtPosEtc, POS_MASK, POS_SHIFT};
use crate::num::ToPrimitive;
use crate::rdbg::STAT_DB;
use derive_more::{Add, Rem, Sub};
use serde::{Deserialize, Serialize};
use std::clone::Clone;
use std::fmt;

// 4 twobits per byte, so unshifted pos is shifted another 2.
const BYTE_SHIFT: u32 = POS_SHIFT + 2;

// twobit shifts are bit positions 0, 2, 4 and 6 in a byte.
const TWOBIT_SHIFT: u32 = POS_SHIFT - 1;

// Two types are defined here.
//
// u64 with only bits set for genomic position, shifted with POS_SHIFT
#[derive(
    Add, Sub, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash, Rem,
)]
pub struct Position(u64);

// only bits set for postion, but not shifted yet
#[derive(Add, Sub, Serialize, Deserialize)]
pub struct BasePos(u64);

/////////////// Position \\\\\\\\\\\\\\\\\\

impl Position {
    pub fn as_u64(&self) -> u64 {
        self.0
    }
    // used both for position0 and 'unset'. The context should make the distinction.
    pub fn zero() -> Self {
        Position(0x0)
    }
    pub fn is_set(&self) -> bool {
        self.0 != 0
    }
    pub fn byte_pos(&self) -> usize {
        // 4 twobits per byte.
        (self.0 >> BYTE_SHIFT) as usize
    }
    // before applying TWOBIT_SHIFT (4 places) I had this:
    // XXX: thread 'main' panicked at 'u32??: TryFromIntError(())', src/new_types/position.rs:49:39
    // TODO: check that distances are sensible.
    pub fn b2_shift(&self) -> u32 {
        let bit_pos = self.0.checked_shr(TWOBIT_SHIFT).unwrap() & 6;
        bit_pos.to_u32().unwrap()
    }

    // for a repetitive mark, return if on period for the distance between current pos (self) and stored
    pub fn get_if_mark_on_period(&self, stored_pos: Position, pd: Position) -> Option<BasePos> {
        let dist = self
            .0
            .checked_sub(stored_pos.0)
            .expect("stored_pos is greater??");
        if dist % pd.0 == 0 {
            Some(BasePos(dist >> POS_SHIFT))
        } else {
            None
        }
    }
    // increment the position one
    pub fn incr(&mut self) {
        self.0 += 1_u64.checked_shl(POS_SHIFT).expect("pos.incr shft");
        dbg_assert_eq!(self.0, self.0 & POS_MASK);
    }
    // increment the position one
    pub fn decr(&mut self) {
        dbg_assert_ne!(self.0, 0);
        self.0 -= 1_u64.checked_shl(POS_SHIFT).expect("pos.incr shft");
    }
}

impl From<BasePos> for Position {
    fn from(base_pos: BasePos) -> Position {
        let shl = u64::from(base_pos).checked_shl(POS_SHIFT).unwrap();
        dbg_assert_eq!(shl, shl & POS_MASK, "basepos bleeds into extension!");
        Position(shl)
    }
}

impl From<ExtPosEtc> for Position {
    fn from(p: ExtPosEtc) -> Position {
        Position(p.as_u64() & POS_MASK)
    }
}

impl fmt::Debug for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:x}", self.0)
    }
}

/////////////// BasePos \\\\\\\\\\\\\\\\\\

impl BasePos {
    pub fn as_u64(&self) -> u64 {
        self.0
    }
    pub fn as_usize(&self) -> usize {
        usize::try_from(self.0).unwrap()
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

impl From<usize> for BasePos {
    fn from(base_pos: usize) -> BasePos {
        BasePos::from(u64::try_from(base_pos).expect("doesn't fit in u64"))
    }
}

impl From<u64> for BasePos {
    fn from(base_pos: u64) -> BasePos {
        let masked = base_pos & (POS_MASK >> POS_SHIFT);
        dbg_assert_eq!(base_pos, masked, "bleeds beyond position!");
        BasePos(base_pos)
    }
}

impl From<BasePos> for u64 {
    fn from(base_pos: BasePos) -> u64 {
        base_pos.0
    }
}
