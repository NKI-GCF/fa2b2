use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use std::fmt;

pub struct ThreeBit(u8);

#[derive(Copy, Clone, PartialEq)]
pub struct TwoBit(u8);

pub struct TwoBitx4(u8);

#[derive(Copy, Clone, PartialEq, Eq)]
pub(super) struct TwoBitDna(u64);

#[derive(Copy, Clone, PartialEq, Eq)]
pub(super) struct TwoBitRcDna(u64);

/// Twobits may be unexpected: N: 0x7, A: 0x0, C: 0x1, T: 0x2, G: 0x3
impl ThreeBit {
    pub fn as_twobit_if_not_n(&self) -> Option<TwoBit> {
        if self.0 < 4 {
            Some(TwoBit(self.0))
        } else {
            None
        }
    }
}

/// No N, 2 bits for code, same as above.
impl TwoBit {
    pub fn pos_shift(&self, pos: Position) -> TwoBitx4 {
        let ob2 = self.0.checked_shl(pos.b2_shift());
        TwoBitx4(ob2.expect("bug shifting from b2"))
    }
    pub fn as_kmer_top(&self, shift: u32) -> u64 {
        (self.0 as u64)
            .checked_shl(shift)
            .expect("bug shifting to u64 top")
    }
    pub fn as_kmer_bottom_rc(&self) -> u64 {
        2 ^ self.0 as u64
    }
}

/// 4 packed twobits per u8.
impl TwoBitx4 {
    pub fn to_b2(&self, pos: Position, for_repeat: bool) -> TwoBit {
        // the third bit (for N) is actually never set, because we don't store those in TwoBitx4
        let ob2 = self.0.checked_shr(pos.b2_shift());
        let b2 = TwoBit(ob2.expect("bug shifting to b2") & 3);
        if !for_repeat {
            dbg_print!("{}: {:?}", BasePos::from(pos).as_u64(), b2);
        }
        b2
    }
    pub fn as_u8(&self) -> u8 {
        self.0
    }
}

impl TwoBitDna {
    /// adds twobit to kmer dna sequences, in the top two bits.
    pub(super) fn add(&mut self, b2: TwoBit, topb2_shift: u32) {
        self.0 = (self.0 >> 2) | b2.as_kmer_top(topb2_shift);
    }
    pub fn as_u64(&self) -> u64 {
        self.0
    }
}

impl TwoBitRcDna {
    /// adds reverse complement of twobit to reverse complement in the bottom.
    pub(super) fn add(&mut self, b2: TwoBit, topb2_shift: u32) {
        self.0 = ((self.0 & ((1 << topb2_shift) - 1)) << 2) | b2.as_kmer_bottom_rc();
    }
    /*fn as_u64(&self) -> u64 {
        self.0
    }*/
}

impl fmt::Debug for TwoBit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            0 => write!(f, "0 (A)"),
            1 => write!(f, "1 (C)"),
            2 => write!(f, "2 (T)"),
            3 => write!(f, "3 (G)"),
            _ => unreachable!(),
        }
    }
}

impl From<&u8> for TwoBitx4 {
    fn from(val: &u8) -> TwoBitx4 {
        TwoBitx4(*val)
    }
}

impl From<&u8> for ThreeBit {
    fn from(val: &u8) -> ThreeBit {
        let b2 = (*val >> 1) & 0x7;
        dbg_print!("{}: {:x}", *val as char, b2);
        ThreeBit(b2)
    }
}
