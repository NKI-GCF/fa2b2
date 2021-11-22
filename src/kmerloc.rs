use crate::rdbg::STAT_DB;
use num::FromPrimitive;
use num_traits::PrimInt;
use std::ops::{AddAssign, SubAssign};

pub trait PriExtPosOri: PrimInt + FromPrimitive + AddAssign + SubAssign {
    fn extension(&self) -> u64;
    fn x(&self) -> usize;
    fn same_ori(&self, p: u64) -> bool;
    fn pos(&self) -> u64;
    fn blacklist(&mut self);
    fn no_pos(&self) -> bool;
    fn extend(&mut self);
    fn set_extension(&mut self, x: u64);
    fn clear_extension(&mut self);
    fn is_same(&self, other: u64) -> bool;
    fn with_ext(&self, x: usize) -> u64;
}

impl PriExtPosOri for u64 {
    fn extension(&self) -> u64 {
        self & 0xFFFF_0000_0000_0000
    }
    fn x(&self) -> usize {
        (self.extension() >> 48) as usize
    }
    fn same_ori(&self, p: u64) -> bool {
        (self & 1) == (p & 1)
    }
    fn pos(&self) -> u64 {
        *self & 0x_FFFF_FFFF_FFFE
    }
    fn blacklist(&mut self) {
        if !self.no_pos() {
            *self &= !0x_FFFF_FFFF_FFFF;
            self.extend();
        }
    }
    fn no_pos(&self) -> bool {
        self.pos() == 0
    }
    fn extend(&mut self) {
        *self += 1 << 48;
    }
    fn set_extension(&mut self, x: u64) {
        dbg_assert!(x <= 0xFFFF);
        *self &= 0xFFFF_FFFF_FFFF;
        *self |= x << 48;
    }
    fn clear_extension(&mut self) {
        *self &= 0x_FFFF_FFFF_FFFF;
    }
    fn is_same(&self, other: u64) -> bool {
        *self == other
    }
    fn with_ext(&self, x: usize) -> u64 {
        self.pos() | ((x as u64) << 48)
    }
}

/// for extended kmers, where position is midpoint
pub trait MidPos: PriExtPosOri {
    fn is_replaceable_by(&self, new_entry: u64) -> bool;
    fn is_set_and_not(&self, other: u64) -> bool;
    fn has_samepos(&self, other: u64) -> bool;
}

impl MidPos for u64 {
    fn has_samepos(&self, other: u64) -> bool {
        self.pos() == other.pos() && {
            // XXX: may want to remove this later if orientation doesn't matter
            dbg_assert!(self.same_ori(other), "{:#x}, {:#x}", self, other);
            true
        }
    }
    fn is_replaceable_by(&self, new_entry: u64) -> bool {
        *self <= new_entry.extension() || self.has_samepos(new_entry)
    }
    fn is_set_and_not(&self, other: u64) -> bool {
        !(self.no_pos() || self.is_same(other))
    }
}

#[derive(new, Clone, PartialEq)]
pub struct KmerLoc<T> {
    pub idx: usize,
    pub p: T,
}
impl<T: PriExtPosOri> KmerLoc<T> {
    pub fn reset(&mut self) {
        self.idx = usize::max_value();
        self.p = T::from_u64(0).unwrap();
    }

    pub fn is_set(&self) -> bool {
        self.idx != usize::max_value()
    }

    pub fn set(&mut self, idx: usize, p: T, x: usize) {
        self.idx = idx;
        self.p = p;
        self.p.set_extension(x as u64);
    }

    /*pub fn priority(&self) -> u64 {
        self.p.to_u64().unwrap() & 0x8000_0000_0000_0000
    }*/
    pub fn next(&mut self, ori: bool, is_template: bool) {
        // ori == true if kmer is for template, then we want 1 in self.p
        let p = self.p.to_u64().unwrap();
        if is_template {
            self.p += T::from_u64(if ori { 3 } else { 2 } - (1 & p)).unwrap();
        } else {
            self.p -= T::from_u64(if ori { 1 } else { 2 } + (1 & p)).unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};
    #[test]
    fn forward() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 100);
        let pick = rng.gen_range(20, 50);
        for _ in 0..pick {
            let ori = rng.gen_range(0, 2);
            kl.next(ori != 0, true);
            dbg_assert_eq!(ori, kl.p & 1);
        }
        dbg_assert_eq!(kl.p & !1, 100 + (pick << 1));
    }
    #[test]
    fn reverse() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 100);
        let pick = rng.gen_range(20, 50);
        for _ in 0..pick {
            let ori = rng.gen_range(0, 2);
            kl.next(ori != 0, false);
            dbg_assert_eq!(ori, kl.p & 1);
        }
        dbg_assert_eq!(kl.p & !1, 100 - (pick << 1));
    }
}
