use crate::kmerloc::ExtPosEtc;
use derive_more::{Add, Rem, Sub};
use serde::{Deserialize, Serialize};
use std::clone::Clone;
use std::convert::TryFrom;
use std::fmt;

const POS_SHIFT: u32 = 4;
const POS_MASK: u64 = 0x00FF_FFFF_FFFF_FFF0;

// 4 twobits per byte, so unshifted pos is shifted another 2.
const BYTE_SHIFT: u32 = POS_SHIFT + 2;

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
    pub fn byte_pos(&self) -> usize {
        // bytepos is calculated before kmer is complete, so we can't assert self.is.set()
        // the strand bit and 2b encoded, so 4 twobits per byte.
        (self.0 >> BYTE_SHIFT) as usize
    }
    //twobit shifts are bit positions 0, 2, 4 and 6 in a byte.
    pub fn b2_shift(&self) -> u32 {
        self.0
            .checked_shr(POS_SHIFT - 1)
            .map(|p| u32::try_from(p).expect("u32??"))
            .unwrap()
            & 6
    }

    // for a repetitive mark, return if on period for the distance between current pos (self) and stored
    pub fn get_if_mark_on_period(&self, stored_pos: Position, pd: u64) -> Option<u64> {
        let dist = self
            .0
            .checked_sub(stored_pos.0)
            .expect("stored_pos is greater??");
        if dist % pd == 0 {
            Some(dist)
        } else {
            None
        }
    }
    // increment the position one
    pub fn incr(&mut self) {
        self.0 += 1_u64.checked_shl(POS_SHIFT).expect("pos.incr shft");
    }
}

impl From<BasePos> for Position {
    fn from(base_pos: BasePos) -> Position {
        Position(u64::from(base_pos).checked_shl(POS_SHIFT).unwrap() & POS_MASK)
    }
}

/*impl From<ExtPosEtc> for Position {
    fn from(p: ExtPosEtc) -> Position {
        Position(p.as_u64() & POS_MASK)
    }
}*/

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

impl From<u64> for BasePos {
    fn from(base_pos: u64) -> BasePos {
        BasePos(base_pos & (POS_MASK >> POS_SHIFT))
    }
}

impl From<BasePos> for u64 {
    fn from(base_pos: BasePos) -> u64 {
        base_pos.0
    }
}

impl From<usize> for BasePos {
    fn from(base_pos: usize) -> BasePos {
        BasePos::from(u64::try_from(base_pos).expect("usize for BasePos doesn't fit in u64"))
    }
}
