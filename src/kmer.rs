use crate::rdbg::STAT_DB;
use num::{FromPrimitive, ToPrimitive, Unsigned};
use num_traits::PrimInt;
use std::cmp;
use std::mem::size_of;
use std::ops::BitXorAssign;

#[derive(Copy, Clone, PartialEq, PartialOrd, Eq)]
/// A kmer that dissociates index and strand orientation
pub struct Kmer<T> {
    pub dna: T,
    pub rc: T,
    topb2: T,
    pub p: u64,
} //^-^\\

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

implement_revcmp!(u8, u16, u32, u64, u128, usize);

impl<T> Kmer<T>
where
    T: Unsigned + FromPrimitive + ToPrimitive + BitXorAssign + PartialOrd + Ord + Eq,
{
    /// get a kmer for this length
    pub fn new(kmerlen: u32, p: u64) -> Self {
        let bitlen = kmerlen * 2;
        let topb2 = bitlen - 2;
        Kmer {
            dna: <T>::zero(),
            rc: <T>::zero(),
            topb2: T::from_u32(topb2).unwrap(),
            p,
        }
    }

    /// adds twobit to kmer sequences, to dna in the top two bits.
    fn add(&mut self, b2: u8) {
        dbg_assert!(b2 <= 3);
        let topb2 = T::to_u64(&self.topb2).unwrap();
        let topless = (1 << topb2) - 1;
        let dna = T::to_u64(&self.dna).unwrap() >> 2;
        let rc = (T::to_u64(&self.rc).unwrap() & topless) << 2;
        self.dna = T::from_u64(dna | (u64::from(b2) << topb2)).unwrap();
        self.rc = T::from_u64(rc ^ 2 ^ u64::from(b2)).unwrap();
        self.p += 2;
    }
    /// true if the kmer is from the template. Palindromes are special.
    pub fn is_template(&self) -> bool {
        self.dna < self.rc || (self.dna == self.rc && (T::to_u64(&self.dna).unwrap() & 1) != 0)
    }

    /// Add twobit to k-mers and return orientation bit as first bit for stored
    pub fn update(&mut self, b2: u8) -> u64 {
        // XXX function is hot
        self.add(b2);
        match self.dna.cmp(&self.rc) {
            cmp::Ordering::Greater => 0,
            cmp::Ordering::Less => 1,
            cmp::Ordering::Equal => T::to_u64(&self.dna).unwrap() & 1,
        }
    }

    /// return an index specific per sequence but the same for the other orientation
    pub fn get_idx(&self, use_min: bool) -> usize {
        let seq = T::to_usize(if (self.dna < self.rc) == use_min {
            &self.dna
        } else {
            &self.rc
        })
        .unwrap();

        // flipped if the top bit is set, to reduce size.
        let overbit = 1 << (T::to_u64(&self.topb2).unwrap() + 1);
        if (seq & overbit) == 0 {
            seq
        } else {
            (overbit - 1) & !seq
        }
    }
    /// return an index specific per sequence but the same for the other orientation
    pub fn get_hash(&self, x: usize) -> (usize, u64) {
        let seq = T::to_usize(if self.dna < self.rc {
            &self.dna
        } else {
            &self.rc
        })
        .unwrap();

        // flipped if the top bit is set, to reduce size.
        let overbit = 1 << (T::to_u64(&self.topb2).unwrap() + 1);
        let m = overbit - 1;
        let idx = if (seq & overbit) == 0 { seq } else { m & !seq };
        (
            idx ^ idx.wrapping_shl(x as u32) & (m ^ ((1 << x) - 1)),
            self.p,
        )
    }
}

impl<T> Ord for Kmer<T>
where
    T: Unsigned + FromPrimitive + ToPrimitive + BitXorAssign + PartialOrd + Ord + Eq,
{
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.get_idx(true).cmp(&other.get_idx(true))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    /// return a new kmer for given index, length and orientation.
    fn kmer_from_idx<T: FromPrimitive>(index: usize, kmerlen: u32, ori: bool) -> Kmer<T> {
        let bitlen = kmerlen * 2;
        let topb2 = bitlen - 2;
        let mut dna = usize::to_u64(&index).unwrap();
        let mut rc = dna.revcmp(kmerlen as usize);
        if dna > rc || (dna == rc && (dna & 1) == 0) {
            let overbit = 1 << (topb2 + 1);
            let overmask = overbit | (overbit - 1);
            dna ^= overmask;
            rc ^= overmask;
        }
        let (dna, rc) = if ori { (dna, rc) } else { (rc, dna) };
        Kmer {
            dna: T::from_u64(dna).unwrap(),
            rc: T::from_u64(rc).unwrap(),
            topb2: T::from_u32(topb2).unwrap(),
        }
    }
    #[test]
    fn test_u64() {
        let mut kmer: Kmer<u64> = Kmer::new(32);
        for i in 0..32 {
            kmer.add(i & 3);
        }
        dbg_assert_eq!(kmer.dna, 0xE4E4E4E4E4E4E4E4); // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        dbg_assert_eq!(kmer.rc, 0xB1B1B1B1B1B1B1B1); // xor 0xaaaaaaaa and reverse per 2bit
        dbg_assert_eq!(kmer.get_idx(true), 0x4E4E4E4E4E4E4E4E); // highest bit is set, so flipped.
    }
    #[test]
    fn test_u32() {
        let mut kmer: Kmer<u32> = Kmer::new(16);
        for i in 0..16 {
            kmer.add(i & 3);
        }
        dbg_assert_eq!(kmer.dna, 0xE4E4E4E4);
        dbg_assert_eq!(kmer.rc, 0xB1B1B1B1);
        dbg_assert_eq!(kmer.get_idx(true), 0x4E4E4E4E);
    }
    #[test]
    fn test_u8() {
        let mut kmer: Kmer<u8> = Kmer::new(4);
        for i in 0..4 {
            kmer.add(i & 3);
        }
        dbg_assert_eq!(kmer.dna, 0xE4);
        dbg_assert_eq!(kmer.rc, 0xB1);
        dbg_assert_eq!(kmer.get_idx(true), 0x4E);
    }
    #[test]
    fn unique() {
        let mut seen = vec![false; 256];
        let mut kmer: Kmer<u8> = Kmer::new(4);
        for i in 0..=255 {
            for j in 0..4 {
                kmer.add((i >> (j << 1)) & 3);
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
        let mut kmer: Kmer<u64> = Kmer::new(kmerlen);
        for _ in 0..32 {
            kmer.add(rng.gen_range(0..4));
        }
        dbg_assert_eq!(kmer.dna.revcmp(kmerlen as usize), kmer.rc);
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

        let mut kmer: Kmer<u64> = Kmer::new(kmerlen);
        for i in 0..last {
            kmer.add(rng.gen_range(0..4));
            if i == pick {
                test_dna = kmer.dna;
                test_rc = kmer.rc;
                test_idx = kmer.get_idx(true);
                test_ori = kmer.dna < kmer.rc;
            }
        }
        dbg_assert_ne!(test_idx, 0xffffffffffffffff);
        let kmer2 = kmer_from_idx(test_idx, kmerlen, test_ori);

        dbg_assert_eq!(test_dna, kmer2.dna);
        dbg_assert_eq!(test_rc, kmer2.rc);

        let kmer3 = kmer_from_idx(test_idx, kmerlen, !test_ori);

        dbg_assert_eq!(test_dna, kmer3.rc);
        dbg_assert_eq!(test_rc, kmer3.dna);
    }
    #[test]
    fn extra() {
        let mut kmer: Kmer<u64> = Kmer::new(4);
        for _ in 0..16 {
            kmer.add(1);
        }
        dbg_assert_eq!(kmer.dna, 0x55);
        dbg_assert_eq!(kmer.rc, 0xff);
    }
}
