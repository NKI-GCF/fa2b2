use crate::rdbg::STAT_DB;
use std::clone::Clone;

pub trait PriExtPosOri: Clone {
    fn no_pos() -> Self;
    fn clear(&mut self);
    fn set_ori(&mut self, p: u64);
    fn get_ori(&self) -> u64;
    fn set(&mut self, p: u64);
    fn get(&self) -> u64;
    fn incr(&mut self);
    fn decr(&mut self);
    fn extension(&self) -> u64;
    fn x(&self) -> usize;
    fn same_ori(&self, p: u64) -> bool;
    fn pos(&self) -> u64;
    fn byte_pos(&self) -> usize;
    fn blacklist(&mut self);
    fn is_no_pos(&self) -> bool;
    fn extend(&mut self);
    fn set_extension(&mut self, x: u64);
    fn clear_extension(&mut self);
    fn is_same(&self, other: u64) -> bool;
    fn with_ext(&self, x: usize) -> u64;
    fn set_dup(&mut self);
    fn is_dup(&self) -> bool;
}

impl PriExtPosOri for u64 {
    fn no_pos() -> Self {
        0x_7FFF_FFFF_FFFE
    }
    fn clear(&mut self) {
        *self = PriExtPosOri::no_pos();
    }
    fn set_ori(&mut self, p: u64) {
        *self ^= (*self ^ p) & 1
    }
    fn get_ori(&self) -> u64 {
        self & 1
    }
    fn set(&mut self, p: u64) {
        *self = p;
    }
    fn get(&self) -> u64 {
        *self
    }
    fn incr(&mut self) {
        *self += 0x2;
    }
    fn decr(&mut self) {
        *self -= 0x2;
    }
    fn extension(&self) -> u64 {
        self & 0xFFFF_0000_0000_0000
    }
    fn x(&self) -> usize {
        (self.extension() >> 48) as usize
    }
    fn same_ori(&self, p: u64) -> bool {
        self.get_ori() == (p & 1)
    }
    fn pos(&self) -> u64 {
        *self & 0x_7FFF_FFFF_FFFE
    }
    fn byte_pos(&self) -> usize {
        // the strand bit and 2b encoded, so 4 twobits per byte.
        (*self & 0x_7FFF_FFFF_FFF8) as usize >> 3
    }
    fn blacklist(&mut self) {
        if !self.is_no_pos() {
            *self &= !0x_7FFF_FFFF_FFFF;
            self.extend();
        }
    }
    fn is_no_pos(&self) -> bool {
        (*self & 0x_7FFF_FFFF_FFFE) == PriExtPosOri::no_pos()
    }
    fn extend(&mut self) {
        *self += 1 << 48;
    }
    fn set_extension(&mut self, x: u64) {
        dbg_assert!(x <= 0xFFFF);
        *self &= 0x7FFF_FFFF_FFFF;
        *self |= x << 48;
    }
    fn clear_extension(&mut self) {
        *self &= 0x_7FFF_FFFF_FFFF;
    }
    fn is_same(&self, other: u64) -> bool {
        *self == other
    }
    fn with_ext(&self, x: usize) -> u64 {
        self.pos() | ((x as u64) << 48)
    }
    fn set_dup(&mut self) {
        *self |= 0x_8000_0000_0000;
    }
    fn is_dup(&self) -> bool {
        *self & 0x_8000_0000_0000 != 0
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
        // blacklisting for smaller extension is setting only extension bits. a value with this
        // extension can be written.
        self.is_no_pos()
            || *self <= new_entry.extension()
            || (self.extension() == new_entry.extension() && self.has_samepos(new_entry))
    }
    fn is_set_and_not(&self, other: u64) -> bool {
        !(self.is_no_pos() || self.is_same(other))
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
        self.p.clear();
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
    /*pub fn next(&mut self, ori: bool, is_template: bool) {
    self.p.set_ori(if ori {1} else {0});*/
    pub fn next(&mut self, p: u64, is_template: bool) {
        if is_template {
            self.p.incr()
        } else {
            self.p.decr()
        }
        self.p.set_ori(p);
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
            kl.next(ori, true);
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
            kl.next(ori, false);
            dbg_assert_eq!(ori, kl.p & 1);
        }
        dbg_assert_eq!(kl.p & !1, 100 - (pick << 1));
    }
}
