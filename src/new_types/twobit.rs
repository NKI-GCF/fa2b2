// (c) Roel Kluin, 2023, GPL v3

use crate::kmerconst::RevCmp;
use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use bitvec::{field::BitField, order::Lsb0, slice::BitSlice, view::BitView};
use derive_more::{From, Into};
use std::fmt;

//TODO: use crate bitvec?

/// N: 0x7, otherwise bits as in [Twobit]
#[derive(Debug)]
// can't leak crate-private type error, unless pub.
/// Twobit: A: 0x0, C: 0x1, T: 0x2, G: 0x3
#[derive(Copy, Clone, PartialEq, new)]
pub struct TwoBit(u8);

/// Four [Twobits] packed in one u8.
pub(crate) struct TwoBitx4(u8);

/// At most 32 [Twobits] packed in one u64, but sometimes not all 32.
#[derive(new, Copy, Clone, PartialEq, Eq, PartialOrd, From, Into)]
pub(crate) struct TwoBitDna(u64);

// TODO: maybe shifting to the reverse complement to the bottom is less performant?
/// At most 32 [Twobits] packed in one u64 and reverse complemented.
/// by convention the reverse complement of a [TwoBitDna].
#[derive(new, Copy, Clone, PartialEq, Eq, From, Into)]
pub(crate) struct TwoBitRcDna(u64);

/// No N, 2 bits for code, same as above.
impl TwoBit {
    // custom from because error is useless
    pub(crate) fn from_u8(val: u8) -> Option<Self> {
        match (val >> 1) & 0x7 {
            b2 if b2 < 4 => Some(TwoBit(b2)),
            _ => None,
        }
    }
    pub(crate) fn as_u8(&self) -> u8 {
        self.0
    }
    pub(crate) fn as_bits(&self) -> &BitSlice<u8, Lsb0> {
        self.0.view_bits::<Lsb0>().split_at(2).0
    }
    pub(crate) fn as_kmer_top(&self, shift: u32) -> u64 {
        (self.0 as u64) << shift
    }
    pub(crate) fn as_kmer_bottom_rc(&self) -> u64 {
        2 ^ self.0 as u64
    }
}

impl From<&BitSlice<u8, Lsb0>> for TwoBit {
    fn from(val: &BitSlice<u8, Lsb0>) -> TwoBit {
        TwoBit(val[0..2].load::<u8>())
    }
}

/// 4 packed twobits per u8.
impl TwoBitx4 {
    pub(crate) fn to_b2(&self, pos: Position) -> TwoBit {
        // the third bit (for N) is actually never set, because we don't store those in TwoBitx4
        let b2 = TwoBit((self.0 >> pos.b2_shift()) & 3);
        dbg_print!("{}: {}", BasePos::from(pos).as_u64(), b2);
        b2
    }
    pub(crate) fn as_u8(&self) -> u8 {
        self.0
    }
}

impl TwoBitDna {
    /// adds twobit to kmer dna sequences, in the top two bits.
    pub(crate) fn add(&mut self, b2: TwoBit, shift: u32) {
        self.0 = (self.0 >> 2) | b2.as_kmer_top(shift);
    }
    #[inline(always)]
    pub(crate) fn lt_strand(&self, rc: TwoBitRcDna) -> bool {
        self.0 < rc.0 || (self.0 == rc.0 && (self.0 & 1) != 0)
    }
    pub(crate) fn to_usize(self) -> usize {
        usize::try_from(self.0).unwrap()
    }
    fn revcmp(&self, k: usize) -> u64 {
        self.0.revcmp(k)
    }
}

impl TwoBitRcDna {
    /// adds reverse complement of twobit to reverse complement in the bottom.
    pub(crate) fn add(&mut self, b2: TwoBit, mask: u64) {
        self.0 = ((self.0 & mask) << 2) | b2.as_kmer_bottom_rc();
    }
    pub(crate) fn to_usize(self) -> usize {
        usize::try_from(self.0).unwrap()
    }
}

impl fmt::Display for TwoBit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            0 => write!(f, "A (0)"),
            1 => write!(f, "C (1)"),
            2 => write!(f, "T (2)"),
            3 => write!(f, "G (3)"),
            _ => unreachable!(),
        }
    }
}

impl From<u8> for TwoBitx4 {
    fn from(val: u8) -> TwoBitx4 {
        TwoBitx4(val)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::XmerHash;
    use crate::kmerconst::{KmerConst, RevCmp};
    use crate::new_types::xmer::Xmer;
    use rand::{thread_rng, Rng};
    use std::cmp;
    /// return a new kmer for given index, length and orientation.
    fn kmer_from_idx(index: usize, kmerlen: usize, ori: bool) -> Xmer {
        let mut dna = index as u64;
        let mut rc = dna.revcmp(kmerlen);
        if dna > rc || (dna == rc && (dna & 1) == 0) {
            let overbit = 1 << (((kmerlen << 1) - 2) + 1);
            let overmask = overbit | (overbit - 1);
            dna ^= overmask;
            rc ^= overmask;
        }
        let (dna, rc) = if ori { (dna, rc) } else { (rc, dna) };
        let mut xmer = Xmer::new();
        xmer.dna = TwoBitDna(dna);
        xmer.rc = TwoBitRcDna(rc);
        xmer
    }

    #[test]
    fn test_u64() {
        let kc = KmerConst::from_bitlen(64, 32, 0);
        let mut xmer: Xmer = Xmer::new();
        for i in 0..32 {
            xmer.update(&kc, TwoBit::new(i & 3));
        }
        assert_eq!(xmer.dna.to_usize(), 0xE4E4E4E4E4E4E4E4); // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        assert_eq!(xmer.rc.to_usize(), 0xB1B1B1B1B1B1B1B1); // xor 0xaaaaaaaa and reverse per 2bit
        let idx = kc.compress_xmer(xmer.get_base_seq());
        assert_eq!(idx, 0x4E4E4E4E4E4E4E4E);
    }
    #[test]
    fn test_u32() {
        let kc = KmerConst::from_bitlen(32, 16, 0);
        let mut xmer: Xmer = Xmer::new();
        for i in 0..16 {
            xmer.update(&kc, TwoBit::new(i & 3));
        }
        assert_eq!(xmer.dna.to_usize(), 0xE4E4E4E4);
        assert_eq!(xmer.rc.to_usize(), 0xB1B1B1B1);
        let idx = kc.compress_xmer(xmer.get_base_seq());
        assert_eq!(idx, 0x4E4E4E4E);
    }
    #[test]
    fn test_u8() {
        let kc = KmerConst::from_bitlen(8, 4, 0);
        let mut xmer: Xmer = Xmer::new();
        for i in 0..4 {
            xmer.update(&kc, TwoBit::new(i & 3));
        }
        assert_eq!(xmer.dna.to_usize(), 0xE4);
        assert_eq!(xmer.rc.to_usize(), 0xB1);
        let idx = kc.compress_xmer(xmer.get_base_seq());
        assert_eq!(idx, 0x4E);
    }
    #[test]
    fn test_revcmp() {
        let mut rng = thread_rng();
        let kmerlen = rng.gen_range(2..32);
        let kc = KmerConst::from_bitlen(kmerlen * 2, kmerlen.into(), 0);
        let kmerlen = usize::from(kmerlen);
        let mut xmer: Xmer = Xmer::new();
        for _ in 0..32 {
            xmer.update(&kc, TwoBit::new(rng.gen_range(0..4)));
        }
        assert_eq!(xmer.dna.revcmp(kmerlen), xmer.rc.0);
    }
    #[test]
    fn test_from_idx() {
        let mut rng = thread_rng();
        let kmerlen = rng.gen_range(12..32);
        let kc = KmerConst::from_bitlen(kmerlen * 2, kmerlen.into(), 0);
        let kmerlen = kmerlen;
        let mut test_dna = 0;
        let mut test_rc = 0;
        let mut test_ori = false;
        let mut test_idx = 0xffffffffffffffff;
        let last = rng.gen_range((kmerlen + 1)..102);
        let pick = rng.gen_range(kmerlen..cmp::max(last - 1, kmerlen + 1));

        let mut xmer: Xmer = Xmer::new();
        for i in 0..last {
            xmer.update(&kc, TwoBit::new(rng.gen_range(0..4)));
            if i == pick {
                test_dna = xmer.dna.to_usize();
                test_rc = xmer.rc.to_usize();
                test_idx = kc.compress_xmer(xmer.get_base_seq()).into();
                test_ori = xmer.dna.lt_strand(xmer.rc);
            }
        }
        assert_ne!(test_idx, 0xffffffffffffffff);
        let xmer2 = kmer_from_idx(test_idx, kmerlen.into(), test_ori);

        assert_eq!(test_dna, xmer2.dna.to_usize());
        assert_eq!(test_rc, xmer2.rc.to_usize());

        let xmer3 = kmer_from_idx(test_idx, kmerlen.into(), !test_ori);

        assert_eq!(test_dna, xmer3.rc.to_usize());
        assert_eq!(test_rc, xmer3.dna.to_usize());
    }
    #[test]
    fn extra() {
        let kc = KmerConst::from_bitlen(8, 4, 0);
        let mut xmer: Xmer = Xmer::new();
        for _ in 0..16 {
            xmer.update(&kc, TwoBit::new(1));
        }
        assert_eq!(xmer.dna.to_usize(), 0x55);
        assert_eq!(xmer.rc.to_usize(), 0xff);
    }
}
