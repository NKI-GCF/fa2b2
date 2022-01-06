use crate::rdbg::STAT_DB;
use std::clone::Clone;
use std::cmp;

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
    fn rep_dup_masked(&self) -> u64;
    fn is_replaceable_by(&self, new_entry: u64) -> bool;
    fn is_set_and_not(&self, other: u64) -> bool;
    fn same_pos_and_ext(&self, new_entry: u64) -> bool;
    fn has_samepos(&self, other: u64) -> bool;
    //fn ext_max() -> Self;
}

impl PriExtPosOri for u64 {
    fn set(&mut self, p: u64) {
        *self = p;
    }
    fn get(&self) -> u64 {
        dbg_assert!(self.is_set());
        *self
    }
    fn pos(&self) -> u64 {
        //assert!(self.is_set());
        *self & 0xF_FFFF_FFFF_FFFE
    }
    fn no_pos() -> Self {
        0x0
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
        dbg_assert!(self.is_set());
        *self ^= (*self ^ p) & 1
    }
    fn get_ori(&self) -> u64 {
        dbg_assert!(self.is_set());
        self & 1
    }
    fn incr(&mut self) {
        dbg_assert!(self.is_set());
        *self += 0x2;
    }
    fn decr(&mut self) {
        dbg_assert!(self.is_set());
        *self -= 0x2;
    }
    // not all extensions may apply, it's dependent on genome size.
    fn extension(&self) -> u64 {
        //assert!(self.is_set());
        self & 0xFF00_0000_0000_0000
    }
    fn x(&self) -> usize {
        (self.extension() >> 56) as usize
    }
    fn same_ori(&self, p: u64) -> bool {
        dbg_assert!(self.is_set());
        self.get_ori() == (p & 1)
    }
    fn byte_pos(&self) -> usize {
        // bytepos is calculated before kmer is complete, so we can't assert self.is.set()
        // the strand bit and 2b encoded, so 4 twobits per byte.
        ((*self & 0xF_FFFF_FFFF_FFFE) >> 3) as usize
    }
    fn blacklist(&mut self) {
        if self.is_set() {
            *self &= !0xF_FFFF_FFFF_FFFF;
            self.extend();
        }
    }
    fn extend(&mut self) {
        dbg_assert!(self.is_set());
        *self += 1 << 56;
    }
    fn set_extension(&mut self, x: u64) {
        dbg_assert!(self.is_set());
        dbg_assert!(x <= 0xFF);
        *self &= 0xF_FFFF_FFFF_FFFF;
        *self |= x << 56;
    }
    fn clear_extension(&mut self) {
        dbg_assert!(self.is_set());
        *self &= 0xF_FFFF_FFFF_FFFF;
    }
    fn is_same(&self, other: u64) -> bool {
        dbg_assert!(self.is_set());
        *self == other
    }
    fn with_ext(&self, x: usize) -> u64 {
        self.pos() | ((x as u64) << 56)
    }
    fn set_dup(&mut self) {
        dbg_assert!(self.is_set());
        *self |= 0x80_0000_0000_0000;
    }
    fn is_dup(&self) -> bool {
        dbg_assert!(self.is_set());
        *self & 0x80_0000_0000_0000 != 0
    }
    fn set_repetitive(&mut self) {
        dbg_assert!(self.is_set());
        *self |= 0x40_0000_0000_0000;
    }
    fn is_repetitive(&self) -> bool {
        dbg_assert!(self.is_set());
        *self & 0x40_0000_0000_0000 != 0
    }
    fn rep_dup_masked(&self) -> u64 {
        dbg_assert!(self.is_set());
        *self & !0xc0_0000_0000_0000
    }
    fn is_replaceable_by(&self, new_entry: u64) -> bool {
        // only extension bits means blacklisting, except for extension 0. pos is always > kmerlen
        new_entry.extension() > self.rep_dup_masked() // TODO: count down extension?
            || (new_entry.extension() == self.extension() && new_entry.pos() <= self.pos())
    }
    fn is_set_and_not(&self, other: u64) -> bool {
        self.is_set() && !self.is_same(other)
    }
    fn same_pos_and_ext(&self, new_entry: u64) -> bool {
        dbg_assert!(self.is_set());
        (*self ^ new_entry) & 0xFF0F_FFFF_FFFF_FFFE == 0
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

#[derive(new, Clone, PartialEq, PartialOrd, Eq)]
pub struct KmerLoc<T>
where
    T: PriExtPosOri + Copy + Ord,
{
    idx: usize,
    pub p: T,
}
impl<T: PriExtPosOri + Copy + Ord> KmerLoc<T> {
    pub fn get(&self) -> (usize, T) {
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

    pub fn set(&mut self, idx: usize, p: T, x: usize) {
        // during rebuilding the strange case occurs that mark is not set, but p is (extension)
        let self_p_extension = self.p.extension();
        let p_extension = p.extension();
        dbg_assert!(self.is_set() || self.p.is_no_pos() || self_p_extension == p_extension);
        self.idx = idx;
        self.p = p;
        self.p.set_extension(x as u64);
    }
}

impl<T> Ord for KmerLoc<T>
where
    T: PriExtPosOri + Copy + Ord,
{
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

    impl<T: PriExtPosOri + Copy> KmerLoc<T> {
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
