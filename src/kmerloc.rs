use crate::rdbg::STAT_DB;
use std::clone::Clone;

pub trait PriExtPosOri: Clone {
    fn set(&mut self, p: u64);
    fn get(&self) -> u64;
    fn pos(&self) -> u64;
    fn no_pos() -> Self;
    fn clear(&mut self);
    fn is_set(&self) -> bool;
    fn is_no_pos(&self) -> bool;
    fn set_ori(&mut self, p: u64);
    fn get_ori(&self) -> u64;
    fn incr(&mut self);
    fn decr(&mut self);
    fn extension(&self) -> u64;
    fn x(&self) -> usize;
    fn same_ori(&self, p: u64) -> bool;
    fn byte_pos(&self) -> usize;
    fn blacklist(&mut self);
    fn extend(&mut self);
    fn set_extension(&mut self, x: u64);
    fn clear_extension(&mut self);
    fn is_same(&self, other: u64) -> bool;
    fn with_ext(&self, x: usize) -> u64;
    fn set_dup(&mut self);
    fn is_dup(&self) -> bool;
    fn set_repetitive(&mut self);
    fn is_repetitive(&self) -> bool;
    //fn ext_max() -> Self;
}

impl PriExtPosOri for u64 {
    fn set(&mut self, p: u64) {
        *self = p;
    }
    fn get(&self) -> u64 {
        *self
    }
    fn pos(&self) -> u64 {
        *self & 0xF_FFFF_FFFF_FFFE
    }
    fn no_pos() -> Self {
        0xF_FFFF_FFFF_FFFE
    }
    fn clear(&mut self) {
        *self = PriExtPosOri::no_pos();
    }
    fn is_set(&self) -> bool {
        self.pos() != PriExtPosOri::no_pos()
    }
    fn is_no_pos(&self) -> bool {
        self.pos() == PriExtPosOri::no_pos()
    }
    fn set_ori(&mut self, p: u64) {
        *self ^= (*self ^ p) & 1
    }
    fn get_ori(&self) -> u64 {
        self & 1
    }
    fn incr(&mut self) {
        *self += 0x2;
    }
    fn decr(&mut self) {
        *self -= 0x2;
    }
    // not all extensions may apply, it's dependent on genome size.
    fn extension(&self) -> u64 {
        self & 0xFF00_0000_0000_0000
    }
    fn x(&self) -> usize {
        (self.extension() >> 56) as usize
    }
    fn same_ori(&self, p: u64) -> bool {
        self.get_ori() == (p & 1)
    }
    fn byte_pos(&self) -> usize {
        // the strand bit and 2b encoded, so 4 twobits per byte.
        self.pos() as usize >> 3
    }
    fn blacklist(&mut self) {
        if self.is_set() {
            *self &= !0xF_FFFF_FFFF_FFFF;
            self.extend();
        }
    }
    fn extend(&mut self) {
        *self += 1 << 56;
    }
    fn set_extension(&mut self, x: u64) {
        dbg_assert!(x <= 0xFF);
        *self &= 0xF_FFFF_FFFF_FFFF;
        *self |= x << 56;
    }
    fn clear_extension(&mut self) {
        *self &= 0xF_FFFF_FFFF_FFFF;
    }
    fn is_same(&self, other: u64) -> bool {
        *self == other
    }
    fn with_ext(&self, x: usize) -> u64 {
        self.pos() | ((x as u64) << 56)
    }
    fn set_dup(&mut self) {
        *self |= 0x80_0000_0000_0000;
    }
    fn is_dup(&self) -> bool {
        *self & 0x80_0000_0000_0000 != 0
    }
    fn set_repetitive(&mut self) {
        *self |= 0x40_0000_0000_0000;
    }
    fn is_repetitive(&self) -> bool {
        *self & 0x40_0000_0000_0000 != 0
    }
    /* does not work because extension is currently limited by readlen.
     * could maybe instead also shift base kmer or bits
     * fn ext_max() -> Self {
        0xFF
    }*/
}

/// for extended kmers, where position is midpoint
pub trait MidPos: PriExtPosOri {
    fn is_replaceable_by(&self, new_entry: u64) -> bool;
    fn is_set_and_not(&self, other: u64) -> bool;
    fn has_samepos(&self, other: u64) -> bool;
    fn same_pos_and_ext(&self, new_entry: u64) -> bool;
}

impl MidPos for u64 {
    fn same_pos_and_ext(&self, new_entry: u64) -> bool {
        (*self ^ new_entry) & 0xFF0F_FFFF_FFFF_FFFE == 0
    }
    fn has_samepos(&self, other: u64) -> bool {
        self.pos() == other.pos() && {
            // XXX: may want to remove this later if orientation doesn't matter
            dbg_assert!(self.same_ori(other), "{:#x}, {:#x}", self, other);
            true
        }
    }
    fn is_replaceable_by(&self, new_entry: u64) -> bool {
        // only extension bits means blacklisting, except for extension 0. pos is always > kmerlen
        self.is_no_pos() || *self <= new_entry.extension() || self.same_pos_and_ext(new_entry)
    }
    fn is_set_and_not(&self, other: u64) -> bool {
        self.is_set() && !self.is_same(other)
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    impl<T: PriExtPosOri> KmerLoc<T> {
        fn next(&mut self, p: u64, is_template: bool) {
            if is_template {
                self.p.incr()
            } else {
                self.p.decr()
            }
            self.p.set_ori(p);
        }
    }
    #[test]
    fn forward() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 100);
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, true);
            dbg_assert_eq!(ori, kl.p & 1);
        }
        dbg_assert_eq!(kl.p & !1, 100 + (pick << 1));
    }
    #[test]
    fn reverse() {
        let mut rng = thread_rng();
        let mut kl = KmerLoc::new(10, 100);
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, false);
            dbg_assert_eq!(ori, kl.p & 1);
        }
        dbg_assert_eq!(kl.p & !1, 100 - (pick << 1));
    }
}
