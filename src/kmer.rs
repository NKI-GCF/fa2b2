use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use num::{FromPrimitive, Unsigned};
use num_traits::PrimInt;
use std::mem::size_of;
use std::{cmp, fmt};

pub struct ThreeBit(u8);

#[derive(PartialEq)]
pub struct TwoBit(u8);

pub struct TwoBitx4(u8);

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
pub struct Kmer<T> {
    pub dna: T,
    pub rc: T,
    pub pos: Position,
    topb2_shift: u32,
} //^-^\\

impl<T> Kmer<T>
where
    T: Unsigned + PrimInt + FromPrimitive,
{
    /// get a kmer for this length
    pub fn new(kmerlen: u32) -> Self {
        let bitlen = kmerlen * 2;
        let topb2_shift = bitlen - 2;
        Kmer {
            dna: <T>::zero(),
            rc: <T>::zero(),
            pos: Position::zero(),
            topb2_shift,
        }
    }

    /// adds twobit to kmer sequences, to dna in the top two bits.
    fn add(&mut self, b2: TwoBit) {
        let kmer_mask = (1 << self.topb2_shift) - 1;
        let dna = T::to_u64(&self.dna).unwrap() >> 2;
        let rc = (T::to_u64(&self.rc).unwrap() & kmer_mask) << 2;
        self.dna = T::from_u64(dna | b2.as_kmer_top(self.topb2_shift)).unwrap();
        self.rc = T::from_u64(rc | b2.as_kmer_bottom_rc()).unwrap();
        self.pos.incr();
    }
    /// true if the kmer is from the template. Palindromes are special.
    pub fn is_template(&self) -> bool {
        self.dna < self.rc || (self.dna == self.rc && (T::to_u64(&self.dna).unwrap() & 1) != 0)
    }

    /// Add twobit to k-mers and return orientation bit as first bit for stored
    pub fn update(&mut self, b2: TwoBit) -> bool {
        // XXX function is hot
        self.add(b2);
        match self.dna.cmp(&self.rc) {
            cmp::Ordering::Greater => false,
            cmp::Ordering::Less => true,
            cmp::Ordering::Equal => T::to_u64(&self.dna).unwrap() & 1 != 0,
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
        let overbit = 1_usize
            .checked_shl(self.topb2_shift + 1)
            .expect("get_idx shft");
        if (seq & overbit) == 0 {
            seq
        } else {
            (overbit - 1) & !seq
        }
    }

    // TODO: implement this instead.
    /// return an index specific per sequence but the same for the other orientation
    pub fn get_hash(&self, x: usize) -> (usize, Position) {
        let seq = T::to_usize(if self.dna < self.rc {
            &self.dna
        } else {
            &self.rc
        })
        .unwrap();

        // flipped if the top bit is set, to reduce size.
        let overbit = 1_usize
            .checked_shl(self.topb2_shift + 1)
            .expect("get_hash shft");
        let m = overbit - 1;
        let idx = if (seq & overbit) == 0 { seq } else { m & !seq };
        (
            idx ^ idx.wrapping_shl(x as u32) & (m ^ ((1 << x) - 1)),
            self.pos,
        )
    }
}

impl<T> PartialOrd for Kmer<T>
where
    T: Unsigned + PrimInt + FromPrimitive,
{
    fn partial_cmp(&self, other: &Kmer<T>) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Ord for Kmer<T>
where
    T: Unsigned + PrimInt + FromPrimitive,
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
    fn kmer_from_idx<T>(index: usize, kmerlen: u32, ori: bool) -> Kmer<T>
    where
        T: Unsigned + PrimInt + FromPrimitive,
    {
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
        kmer.dna = T::from_u64(dna).unwrap();
        kmer.rc = T::from_u64(rc).unwrap();
        kmer
    }
    #[test]
    fn test_u64() {
        let mut kmer: Kmer<u64> = Kmer::new(32);
        for i in 0..32 {
            kmer.add(TwoBit(i & 3));
        }
        dbg_assert_eq!(kmer.dna, 0xE4E4E4E4E4E4E4E4); // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        dbg_assert_eq!(kmer.rc, 0xB1B1B1B1B1B1B1B1); // xor 0xaaaaaaaa and reverse per 2bit
        dbg_assert_eq!(kmer.get_idx(true), 0x4E4E4E4E4E4E4E4E); // highest bit is set, so flipped.
    }
    #[test]
    fn test_u32() {
        let mut kmer: Kmer<u32> = Kmer::new(16);
        for i in 0..16 {
            kmer.add(TwoBit(i & 3));
        }
        dbg_assert_eq!(kmer.dna, 0xE4E4E4E4);
        dbg_assert_eq!(kmer.rc, 0xB1B1B1B1);
        dbg_assert_eq!(kmer.get_idx(true), 0x4E4E4E4E);
    }
    #[test]
    fn test_u8() {
        let mut kmer: Kmer<u8> = Kmer::new(4);
        for i in 0..4 {
            kmer.add(TwoBit(i & 3));
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
        let mut kmer: Kmer<u64> = Kmer::new(kmerlen);
        for _ in 0..32 {
            kmer.add(TwoBit(rng.gen_range(0..4)));
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
            kmer.add(TwoBit(rng.gen_range(0..4)));
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
            kmer.add(TwoBit(1));
        }
        dbg_assert_eq!(kmer.dna, 0x55);
        dbg_assert_eq!(kmer.rc, 0xff);
    }
}
