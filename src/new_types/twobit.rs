use crate::kmerconst::RevCmp;
use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use derive_more::{From, Into};
use std::fmt;

pub struct ThreeBit(u8);

#[derive(Copy, Clone, PartialEq, new)]
pub struct TwoBit(u8);

pub struct TwoBitx4(u8);

#[derive(new, Copy, Clone, PartialEq, Eq, PartialOrd, From, Into)]
pub(super) struct TwoBitDna(u64);

#[derive(new, Copy, Clone, PartialEq, Eq, From, Into)]
pub(super) struct TwoBitRcDna(u64);

/// Twobits may be unexpected: N: 0x7, A: 0x0, C: 0x1, T: 0x2, G: 0x3
impl ThreeBit {
    pub(crate) fn as_twobit_if_not_n(&self) -> Option<TwoBit> {
        if self.0 < 4 {
            Some(TwoBit(self.0))
        } else {
            None
        }
    }
}

/// No N, 2 bits for code, same as above.
impl TwoBit {
    pub(crate) fn pos_shift(&self, pos: Position) -> TwoBitx4 {
        let ob2 = self.0.checked_shl(pos.b2_shift());
        TwoBitx4(ob2.expect("bug shifting from b2"))
    }
    pub(crate) fn as_kmer_top(&self, shift: u32) -> u64 {
        (self.0 as u64)
            .checked_shl(shift)
            .expect("bug shifting to u64 top")
    }
    pub(crate) fn as_kmer_bottom_rc(&self) -> u64 {
        2 ^ self.0 as u64
    }
}

/// 4 packed twobits per u8.
impl TwoBitx4 {
    pub(crate) fn to_b2(&self, pos: Position, for_repeat: bool) -> TwoBit {
        // the third bit (for N) is actually never set, because we don't store those in TwoBitx4
        let ob2 = self.0.checked_shr(pos.b2_shift());
        let b2 = TwoBit(ob2.expect("bug shifting to b2") & 3);
        if !for_repeat {
            dbg_print!("{}: {:?}", BasePos::from(pos).as_u64(), b2);
        }
        b2
    }
    pub(crate) fn as_u8(&self) -> u8 {
        self.0
    }
}

impl TwoBitDna {
    /// adds twobit to kmer dna sequences, in the top two bits.
    pub(super) fn add(&mut self, b2: TwoBit, topb2_shift: u32) {
        self.0 = (self.0 >> 2) | b2.as_kmer_top(topb2_shift);
    }
    /*pub(crate) fn as_u64(&self) -> u64 {
        self.0
    }*/
    pub(crate) fn lt_strand(&self, rc: TwoBitRcDna) -> bool {
        self.0 < rc.0 || (self.0 == rc.0 && (self.0 & 1) != 0)
    }
    pub(crate) fn to_usize(self) -> usize {
        usize::try_from(self.0).unwrap()
    }
    fn revcmp(&self, k: u32) -> u64 {
        self.0.revcmp(k)
    }
}

impl TwoBitRcDna {
    /// adds reverse complement of twobit to reverse complement in the bottom.
    pub(super) fn add(&mut self, b2: TwoBit, topb2_shift: u32) {
        self.0 = ((self.0 & ((1 << topb2_shift) - 1)) << 2) | b2.as_kmer_bottom_rc();
    }
    pub(crate) fn to_usize(self) -> usize {
        usize::try_from(self.0).unwrap()
    }

    /*fn as_u64(&self) -> u64 {
        self.0
    }*/
}

impl fmt::Debug for TwoBit {
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

impl From<&u8> for TwoBitx4 {
    fn from(val: &u8) -> TwoBitx4 {
        TwoBitx4(*val)
    }
}

impl From<(Position, u8)> for ThreeBit {
    fn from(val: (Position, u8)) -> ThreeBit {
        let b2 = (val.1 >> 1) & 0x7;
        dbg_print!(
            "{}: {} ({:x})",
            BasePos::from(val.0).as_u64(),
            val.1 as char,
            b2
        );
        ThreeBit(b2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::{KmerConst, RevCmp};
    use crate::new_types::xmer::xmer_hash;
    use crate::new_types::xmer::Xmer;
    use rand::{thread_rng, Rng};
    use std::cmp;

    /// return a new kmer for given index, length and orientation.
    fn kmer_from_idx(index: usize, kmerlen: u32, ori: bool) -> Xmer {
        let mut dna = index as u64;
        let mut rc = dna.revcmp(kmerlen);
        if dna > rc || (dna == rc && (dna & 1) == 0) {
            let overbit = 1 << (((kmerlen << 1) - 2) + 1);
            let overmask = overbit | (overbit - 1);
            dna ^= overmask;
            rc ^= overmask;
        }
        let (dna, rc) = if ori { (dna, rc) } else { (rc, dna) };
        let mut kmer = Xmer::new(kmerlen);
        kmer.dna = TwoBitDna(dna);
        kmer.rc = TwoBitRcDna(rc);
        kmer
    }
    #[test]
    fn xmer_no_flip() {
        let mut kmer: Xmer = Xmer::new(32);
        assert_eq!(kmer.kmerlen, 32);

        for i in 0..32 {
            kmer.update(TwoBit::new(i & 3));
        }
        let dna = 0xE4E4_E4E4_E4E4_E4E4; // GTCA.. => 3210.. (in 2bits)
        let rc = dna.revcmp(kmer.kmerlen); // xor 0xaaaaaaaa and reverse per 2bit
        assert_eq!(0xB1B1_B1B1_B1B1_B1B1_u64, rc); // TGAC.. 2301.. (in 2bits)

        assert_eq!(kmer.dna.0, dna);
        assert_eq!(kmer.rc.0, rc);

        let mark = kmer.get_hash_and_p(0x55);

        // lowest is rc, so that will be hashed by above function.
        let test = xmer_hash(rc as usize, 0x55, kmer.kmerlen) as u64;
        assert_eq!(test, 0x11B1B1B1A0_u64, "{:#X}", test);

        assert_eq!(test & 0x8000_0000_0000_0000, 0);
        // its highest 'overbit' is not set, so no bitwise complementing occurs.

        assert_eq!(mark.0 as u64, test, "test: {:#x} {:#x}", mark.0, test);

        let top = (rc & !0x55 & 0xFFFF_FFFF) << 32;
        let bottom = (rc >> 32) & 0x55;
        assert_eq!(
            mark.0 as u64,
            rc ^ (top | bottom),
            "top|bottom: {:#x} {:#x}",
            mark.0 as u64,
            rc ^ (top | bottom)
        );
        assert_eq!(mark.0 as u64, 0x11B1B1B1A0);
        assert_eq!(mark.1.x(), 0x55);

        let undo = xmer_hash(mark.0, 0x55, kmer.kmerlen) as u64;
        assert_eq!(undo, rc, "undo: {:#x} {:#x}", undo, rc);
    }

    #[test]
    fn xmer_with_flip() {
        let kc = KmerConst::from_bitlen(64);
        let mut kmer: Xmer = Xmer::new(kc.kmerlen.try_into().unwrap());
        assert_eq!(kmer.kmerlen as usize, kc.kmerlen);
        assert_eq!(kmer.kmerlen, 32);

        let dna = 0x0147_3840_FABC_025F_u64;
        for i in 0..32 {
            kmer.update(TwoBit::new(((dna >> (i << 1)) as u8) & 3));
        }
        let rc = dna.revcmp(kmer.kmerlen); // xor 0xaaaaaaaa and reverse per 2bit
        assert_eq!(0x5F2A_9405_AB86_7BEA, rc, "{:#X}", rc); // TGAC.. 2301.. (in 2bits)

        assert_eq!(kmer.dna.0, dna); // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        assert_eq!(kmer.rc.0, rc); // xor 0xaaaaaaaa and reverse per 2bit

        let mark = kmer.get_hash_and_p(0x55);

        // lowest is dna, so that will be hashed by above function.
        let mut test = xmer_hash(dna as usize, 0x55, kmer.kmerlen) as u64;
        assert_eq!(test, 0xFBFB_3A4A_FABC_021F, "{:#X}", test);

        // its highest 'overbit' is set, so bitwise complementing occurs.
        assert_ne!(test & 0x8000_0000_0000_0000, 0, "{:#x}", test);
        test = !test & 0x7FFF_FFFFFFFF_FFFF;

        assert_eq!(mark.0 as u64, test, "test: {:#x} {:#x}", mark.0, test);

        let top = (dna & !0x55 & 0xFFFF_FFFF) << 32;
        let bottom = (dna >> 32) & 0x55;
        assert_eq!(
            mark.0 as u64,
            (!dna ^ (top | bottom)) & 0x7FFF_FFFFFFFF_FFFF,
            "top|bottom: {:#x} {:#x}",
            mark.0 as u64,
            (!dna ^ (top | bottom)) & 0x7FFF_FFFFFFFF_FFFF
        );
        assert_eq!(mark.0 as u64, 0x404c5b50543fde0, "{:#x}", mark.0 as u64);
        assert_eq!(mark.1.x(), 0x55);

        let mut undo = mark.0 as u64;
        undo ^= 0xFFFF_FFFFFFFF_FFFF;
        undo = xmer_hash(undo as usize, 0x55, kmer.kmerlen) as u64;
        assert_eq!(undo, dna, "undo: {:#x} {:#x}", undo, dna);

        let mut next = kc.get_next_xmer(mark.0, mark.1).unwrap();
        assert_eq!(next.1.x(), 0x56);

        let same = kmer.get_hash_and_p(0x56);
        assert_eq!(next.0, same.0, "undo: {:#x} {:#x}", next.0, same.0);

        for _ in 0..100 {
            next = kc.get_next_xmer(next.0, next.1).unwrap();
            let same = kmer.get_hash_and_p(next.1.x());
            assert_eq!(next.0, same.0, "undo: {:#x} {:#x}", next.0, same.0);
        }
    }

    #[test]
    fn test_u64() {
        let mut kmer: Xmer = Xmer::new(32);
        for i in 0..32 {
            kmer.update(TwoBit::new(i & 3));
        }
        assert_eq!(kmer.dna.0, 0xE4E4E4E4E4E4E4E4); // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        assert_eq!(kmer.rc.0, 0xB1B1B1B1B1B1B1B1); // xor 0xaaaaaaaa and reverse per 2bit
    }
    #[test]
    fn test_u32() {
        let mut kmer: Xmer = Xmer::new(16);
        for i in 0..16 {
            kmer.update(TwoBit::new(i & 3));
        }
        assert_eq!(kmer.dna.0, 0xE4E4E4E4);
        assert_eq!(kmer.rc.0, 0xB1B1B1B1);
        let idx: usize = kmer.get_baseidx().into();
        assert_eq!(idx, 0x4E4E4E4E);
    }
    #[test]
    fn test_u8() {
        let mut kmer: Xmer = Xmer::new(4);
        for i in 0..4 {
            kmer.update(TwoBit::new(i & 3));
        }
        assert_eq!(kmer.dna.0, 0xE4);
        assert_eq!(kmer.rc.0, 0xB1);
        let idx: usize = kmer.get_baseidx().into();
        assert_eq!(idx, 0x4E);
    }
    #[test]
    fn unique() {
        let mut seen = vec![false; 256];
        let mut kmer: Xmer = Xmer::new(4);
        for i in 0..=255 {
            for j in 0..4 {
                kmer.update(TwoBit::new((i >> (j << 1)) & 3));
            }
            let idx: usize = kmer.get_baseidx().into();
            let x = (if kmer.is_template() { 1 } else { 0 }) | idx << 1;
            assert!(!seen[x], "0x{:x} already seen!", x);
            seen[x] = true;
        }
        assert_eq!(vec![true; 256], seen);
    }
    #[test]
    fn test_revcmp() {
        let mut rng = thread_rng();
        let kmerlen = rng.gen_range(2..32);
        let mut kmer: Xmer = Xmer::new(kmerlen);
        for _ in 0..32 {
            kmer.update(TwoBit::new(rng.gen_range(0..4)));
        }
        assert_eq!(kmer.dna.revcmp(kmerlen), kmer.rc.0);
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

        let mut kmer: Xmer = Xmer::new(kmerlen);
        for i in 0..last {
            kmer.update(TwoBit::new(rng.gen_range(0..4)));
            if i == pick {
                test_dna = kmer.dna.0;
                test_rc = kmer.rc.0;
                test_idx = kmer.get_baseidx().into();
                test_ori = kmer.dna.lt_strand(kmer.rc);
            }
        }
        assert_ne!(test_idx, 0xffffffffffffffff);
        let kmer2 = kmer_from_idx(test_idx, kmerlen, test_ori);

        assert_eq!(test_dna, kmer2.dna.0);
        assert_eq!(test_rc, kmer2.rc.0);

        let kmer3 = kmer_from_idx(test_idx, kmerlen, !test_ori);

        assert_eq!(test_dna, kmer3.rc.0);
        assert_eq!(test_rc, kmer3.dna.0);
    }
    #[test]
    fn extra() {
        let mut kmer: Xmer = Xmer::new(4);
        for _ in 0..16 {
            kmer.update(TwoBit::new(1));
        }
        assert_eq!(kmer.dna.0, 0x55);
        assert_eq!(kmer.rc.0, 0xff);
    }
}
