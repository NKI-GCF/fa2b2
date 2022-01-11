use crate::kmerloc::ExtPosEtc;
use crate::new_types::extension::Extension;
use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use std::{cmp, fmt};

pub struct ThreeBit(u8);

#[derive(PartialEq)]
pub struct TwoBit(u8);

pub struct TwoBitx4(u8);

#[derive(Copy, Clone, PartialEq, Eq)]
pub struct TwoBitDna(u64);

#[derive(Copy, Clone, PartialEq, Eq)]
pub struct TwoBitRcDna(u64);

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
    fn add(&mut self, b2: TwoBit, topb2_shift: u32) {
        self.0 = (self.0 >> 2) | b2.as_kmer_top(topb2_shift);
    }
    pub fn as_u64(&self) -> u64 {
        self.0
    }
}

impl TwoBitRcDna {
    /// adds reverse complement of twobit to reverse complement in the bottom.
    fn add(&mut self, b2: TwoBit, topb2_shift: u32) {
        self.0 = ((self.0 & ((1 << topb2_shift) - 1)) << 2) | b2.as_kmer_bottom_rc();
    }
    fn as_u64(&self) -> u64 {
        self.0
    }
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

#[derive(Copy, Clone, PartialEq, Eq)]
/// A kmer that dissociates index and strand orientation
pub struct Kmer {
    pub dna: TwoBitDna,
    pub rc: TwoBitRcDna,
    pub pos: Position,
    kmerlen: u32,
    topb2_shift: u32,
} //^-^\\

impl Kmer {
    /// get a kmer for this length
    pub fn new(kmerlen: u32) -> Self {
        let bitlen = kmerlen * 2;
        let topb2_shift = bitlen - 2;
        Kmer {
            dna: TwoBitDna(0),
            rc: TwoBitRcDna(0),
            pos: Position::zero(),
            kmerlen,
            topb2_shift,
        }
    }

    /// adds twobit to kmer sequences, to dna in the top two bits.
    fn add(&mut self, b2: TwoBit) {
        let kmer_mask = (1 << self.topb2_shift) - 1;
        let dna = self.dna.0 >> 2;
        let rc = (self.rc.0 & kmer_mask) << 2;
        self.dna.0 = dna | b2.as_kmer_top(self.topb2_shift);
        self.rc.0 = rc | b2.as_kmer_bottom_rc();
        self.pos.incr();
    }
    /// true if the kmer is from the template. Palindromes are special.
    pub fn is_template(&self) -> bool {
        self.dna.0 < self.rc.0 || (self.dna.0 == self.rc.0 && (self.dna.0 & 1) != 0)
    }

    /// Add twobit to k-mers and return orientation bit as first bit for stored
    pub fn update(&mut self, b2: TwoBit) -> bool {
        // XXX function is hot
        self.add(b2);
        match self.dna.0.cmp(&self.rc.0) {
            cmp::Ordering::Greater => false,
            cmp::Ordering::Less => true,
            cmp::Ordering::Equal => self.dna.0 & 1 != 0,
        }
    }
    //TODO: remove use_min
    pub fn get_base_seq(&self, use_min: bool) -> u64 {
        if (self.dna.0 < self.rc.0) == use_min {
            self.dna.0
        } else {
            self.rc.0
        }
    }

    /// return an index specific per sequence but the same for the other orientation
    pub fn get_idx(&self, use_min: bool) -> usize {
        let seq = usize::try_from(self.get_base_seq(use_min)).unwrap();

        // flipped if the top bit is set, to reduce size.
        let overbit = 1 << (self.topb2_shift + 1);
        if (seq & overbit) == 0 {
            seq
        } else {
            (overbit - 1) & !seq
        }
    }

    // TODO: implement this instead.
    /// return an index specific per sequence but the same for the other orientation
    pub fn get_hash_and_p(&self, x: usize) -> (usize, ExtPosEtc) {
        let t = 1 << self.kmerlen;
        dbg_assert!(x < t);

        let mut hash = usize::try_from(self.get_base_seq(true)).unwrap();

        hash ^= ((hash & !x & (t - 1)) << self.kmerlen) | ((hash >> self.kmerlen) & x);

        let overbit = 1_usize << (self.topb2_shift + 1);
        let p = ExtPosEtc::from((Extension::from(x), self.pos));

        if hash & overbit != 0 {
            (hash, p)
        } else {
            ((overbit - 1) & !hash, p)
        }
    }
}

impl PartialOrd for Kmer {
    fn partial_cmp(&self, other: &Kmer) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Kmer {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.get_idx(true).cmp(&other.get_idx(true))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num::{FromPrimitive, PrimInt, Unsigned};
    use rand::{thread_rng, Rng};
    use std::mem::size_of;

    pub trait RevCmp<T: PrimInt + FromPrimitive> {
        fn revcmp(self, kmerlen: usize) -> T;
    }

    /// create bitmask. e.g. dvm::<u32>(0xf0  0xff) => 0xf0_f0_f0_f0
    #[inline]
    fn dvm<T: PrimInt + FromPrimitive>(numerator: u32, divisor: u32) -> T {
        let base = T::max_value() / T::from_u32(divisor).unwrap();
        T::from_u32(numerator).unwrap() * base
    }

    macro_rules! implement_revcmp { ($($ty:ty),*) => ($(
        /// give twobit reverse complent for given kmerlen
        impl RevCmp<$ty> for $ty {
            #[inline]
            fn revcmp(self, kmerlen: usize) -> $ty {
                let mut seq = self.swap_bytes() ^ dvm::<$ty>(2, 3);

                let left_nibbles = (seq & dvm::<$ty>(0xf, 0xff)) << 4;
                let right_nibbles = (seq & dvm::<$ty>(0xf0, 0xff)) >> 4;
                seq = left_nibbles | right_nibbles;

                let all_left_two_bits = (seq & dvm::<$ty>(0x3, 0xf)) << 2;
                let all_right_two_bits = (seq & dvm::<$ty>(0xc, 0xf)) >> 2;
                seq = all_left_two_bits | all_right_two_bits;

                seq >> (size_of::<$ty>() * 8 - kmerlen * 2)
            }
        }
        )*)
    }
    implement_revcmp!(u64);

    /// return a new kmer for given index, length and orientation.
    fn kmer_from_idx(index: usize, kmerlen: u32, ori: bool) -> Kmer {
        let mut dna = index as u64;
        let mut rc = dna.revcmp(kmerlen as usize);
        if dna > rc || (dna == rc && (dna & 1) == 0) {
            let overbit = 1 << (((kmerlen << 1) - 2) + 1);
            let overmask = overbit | (overbit - 1);
            dna ^= overmask;
            rc ^= overmask;
        }
        let (dna, rc) = if ori { (dna, rc) } else { (rc, dna) };
        let mut kmer = Kmer::new(kmerlen);
        kmer.dna.0 = dna;
        kmer.rc.0 = rc;
        kmer
    }
    #[test]
    fn test_u64() {
        let mut kmer: Kmer = Kmer::new(32);
        for i in 0..32 {
            kmer.add(TwoBit(i & 3));
        }
        dbg_assert_eq!(kmer.dna.0, 0xE4E4E4E4E4E4E4E4); // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        dbg_assert_eq!(kmer.rc.0, 0xB1B1B1B1B1B1B1B1); // xor 0xaaaaaaaa and reverse per 2bit
        dbg_assert_eq!(kmer.get_idx(true), 0x4E4E4E4E4E4E4E4E); // highest bit is set, so flipped.
    }
    #[test]
    fn test_u32() {
        let mut kmer: Kmer = Kmer::new(16);
        for i in 0..16 {
            kmer.add(TwoBit(i & 3));
        }
        dbg_assert_eq!(kmer.dna.0, 0xE4E4E4E4);
        dbg_assert_eq!(kmer.rc.0, 0xB1B1B1B1);
        dbg_assert_eq!(kmer.get_idx(true), 0x4E4E4E4E);
    }
    #[test]
    fn test_u8() {
        let mut kmer: Kmer = Kmer::new(4);
        for i in 0..4 {
            kmer.add(TwoBit(i & 3));
        }
        dbg_assert_eq!(kmer.dna.0, 0xE4);
        dbg_assert_eq!(kmer.rc.0, 0xB1);
        dbg_assert_eq!(kmer.get_idx(true), 0x4E);
    }
    #[test]
    fn unique() {
        let mut seen = vec![false; 256];
        let mut kmer: Kmer = Kmer::new(4);
        for i in 0..=255 {
            for j in 0..4 {
                kmer.add(TwoBit((i >> (j << 1)) & 3));
            }
            let x = (if kmer.is_template() { 1 } else { 0 }) | kmer.get_idx(true) << 1;
            dbg_assert!(!seen[x], "0x{:x} already seen!", x);
            seen[x] = true;
        }
        dbg_assert_eq!(vec![true; 256], seen);
    }
    #[test]
    fn test_revcmp() {
        let mut rng = thread_rng();
        let kmerlen = rng.gen_range(2..32);
        let mut kmer: Kmer = Kmer::new(kmerlen);
        for _ in 0..32 {
            kmer.add(TwoBit(rng.gen_range(0..4)));
        }
        dbg_assert_eq!(kmer.dna.0.revcmp(kmerlen as usize), kmer.rc.0);
    }
    #[test]
    fn test_from_idx() {
        let mut rng = thread_rng();
        let kmerlen = rng.gen_range(12..32);
        let mut test_dna = 0;
        let mut test_rc = 0;
        let mut test_ori = false;
        let mut test_idx = 0xffffffffffffffff;
        let last = rng.gen_range((kmerlen + 1)..102);
        let pick = rng.gen_range(kmerlen..cmp::max(last - 1, kmerlen + 1));

        let mut kmer: Kmer = Kmer::new(kmerlen);
        for i in 0..last {
            kmer.add(TwoBit(rng.gen_range(0..4)));
            if i == pick {
                test_dna = kmer.dna.0;
                test_rc = kmer.rc.0;
                test_idx = kmer.get_idx(true);
                test_ori = kmer.dna.0 < kmer.rc.0;
            }
        }
        dbg_assert_ne!(test_idx, 0xffffffffffffffff);
        let kmer2 = kmer_from_idx(test_idx, kmerlen, test_ori);

        dbg_assert_eq!(test_dna, kmer2.dna.0);
        dbg_assert_eq!(test_rc, kmer2.rc.0);

        let kmer3 = kmer_from_idx(test_idx, kmerlen, !test_ori);

        dbg_assert_eq!(test_dna, kmer3.rc.0);
        dbg_assert_eq!(test_rc, kmer3.dna.0);
    }
    #[test]
    fn extra() {
        let mut kmer: Kmer = Kmer::new(4);
        for _ in 0..16 {
            kmer.add(TwoBit(1));
        }
        dbg_assert_eq!(kmer.dna.0, 0x55);
        dbg_assert_eq!(kmer.rc.0, 0xff);
    }
}
