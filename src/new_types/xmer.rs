use super::{
    extension::Extension,
    position::Position,
    twobit::{TwoBit, TwoBitDna, TwoBitRcDna},
};
use crate::kmerconst::KmerConst;
use crate::kmerconst::XmerHash;
use crate::new_types::extended_position::ExtPosEtc;
use crate::xmer_location::XmerLoc;
use std::cmp;
//use crate::rdbg::STAT_DB;

#[derive(Copy, Clone, PartialEq, Eq)]
/// An xmer is a kmer that dissociates strand orientation.
pub struct Xmer {
    pub(super) dna: TwoBitDna,
    pub(super) rc: TwoBitRcDna,
    pub(crate) pos: Position,
    base_index: usize,
} //^-^\\

impl Xmer {
    /// get an xmer with no bits set and position 0
    pub(crate) fn new() -> Self {
        Xmer {
            dna: TwoBitDna::new(0),
            rc: TwoBitRcDna::new(0),
            pos: Position::default(),
            base_index: usize::MAX,
        }
    }

    /// true if the xmer is from the template. Palindromes are special.
    #[inline(always)]
    pub(crate) fn is_template(&self) -> bool {
        self.dna.lt_strand(self.rc)
    }

    /// adds twobit to k-mer sequences, to dna in the top two bits. Returns orientation.
    pub(crate) fn update(&mut self, kc: &KmerConst, b2: TwoBit) -> bool {
        self.dna.add(b2, kc.dna_topb2_shift);
        self.rc.add(b2, kc.rc_mask);
        self.pos.incr();
        // XXX: why not make pos a ExtPosEtc and store orientation is first bit?
        let orientation = self.is_template();
        let seq = if orientation {
            self.dna.to_usize()
        } else {
            self.rc.to_usize()
        };

        self.base_index = kc.xmer_hash(seq, kc.seed);

        if self.base_index & kc.overbit != 0 {
            self.base_index = kc.xmer_mask & !self.base_index;
        }
        orientation
    }
    #[inline(always)]
    pub(super) fn get_base_seq(&self) -> usize {
        if self.is_template() {
            self.dna.to_usize()
        } else {
            self.rc.to_usize()
        }
    }
    /// an extension specific index
    pub(crate) fn get_hash(&self, kc: &KmerConst, x: usize) -> usize {
        kc.hash_and_compress(self.get_base_seq(), x)
    }
    pub(crate) fn get_hash_and_p(&self, kc: &KmerConst, x: usize) -> XmerLoc {
        let p = ExtPosEtc::from((Extension::from(x), self.pos));
        XmerLoc::new(self.get_hash(kc, x), p)
    }
}

impl PartialOrd for Xmer {
    fn partial_cmp(&self, other: &Xmer) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Xmer {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.base_index.cmp(&other.base_index)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmerconst::RevCmp;
    use crate::new_types::twobit::TwoBit;
    use crate::new_types::xmer::Xmer;

    #[test]
    fn unique() {
        let kc = KmerConst::from_bitlen(8, 4, 0);
        let mut seen = vec![false; 256];
        let mut xmer: Xmer = Xmer::new();
        for i in 0..=255 {
            for j in 0..4 {
                xmer.update(&kc, TwoBit::new((i >> (j << 1)) & 3));
            }
            let idx = kc.compress_xmer(xmer.get_base_seq());
            let x = (if xmer.is_template() { 1 } else { 0 }) | idx << 1;
            assert!(!seen[x], "0x{:x} already seen!", x);
            seen[x] = true;
        }
        assert_eq!(vec![true; 256], seen);
    }
    #[test]
    fn xmer_reversability_in_xmer_wo_flip() {
        let kc = KmerConst::from_bitlen(64, 32, 0);
        let mut xmer: Xmer = Xmer::new();
        assert_eq!(kc.kmerlen, 32);

        for i in 0..32 {
            xmer.update(&kc, TwoBit::new(i & 3));
        }
        let dna = 0xE4E4_E4E4_E4E4_E4E4; // GTCA.. => 3210.. (in 2bits)
        let rc = dna.revcmp(kc.kmerlen.try_into().unwrap()); // xor 0xaaaaaaaa and reverse per 2bit
        assert_eq!(0xB1B1_B1B1_B1B1_B1B1, rc); // TGAC.. 2301.. (in 2bits)

        assert_eq!(xmer.dna.to_usize(), dna);
        assert_eq!(xmer.rc.to_usize(), rc);

        let mark = xmer.get_hash_and_p(&kc, 0x55);

        // lowest is rc, so that will be hashed by above function.
        let test = kc.xmer_hash(rc, 0x55);
        assert_eq!(test, 0x11B1B1B1A0, "{:#X}", test);

        assert_eq!(test & 0x8000_0000_0000_0000, 0);
        // its highest 'overbit' is not set, so no bitwise complementing occurs.

        assert_eq!(mark.idx, test, "test: {:#x} {:#x}", mark.idx, test);

        let top = (rc & !0x55 & 0xFFFF_FFFF) << 32;
        let bottom = (rc >> 32) & 0x55;
        assert_eq!(
            mark.idx,
            rc ^ (top | bottom),
            "top|bottom: {:#x} {:#x}",
            mark.idx,
            rc ^ (top | bottom)
        );
        assert_eq!(mark.idx, 0x11B1B1B1A0);
        assert_eq!(mark.p.x(), 0x55);

        let undo = kc.xmer_hash(mark.idx, 0x55);
        assert_eq!(undo, rc, "undo: {:#x} {:#x}", undo, rc);
    }

    #[test]
    fn xmer_reversability_in_xmer_with_flip() {
        let kc = KmerConst::from_bitlen(64, 32, 0);
        let mut xmer: Xmer = Xmer::new();
        assert_eq!(kc.kmerlen, kc.kmerlen);
        assert_eq!(kc.kmerlen, 32);

        let dna = 0x0147_3840_FABC_025F;
        for i in 0..32 {
            xmer.update(&kc, TwoBit::new(((dna >> (i << 1)) as u8) & 3));
        }
        let rc = dna.revcmp(kc.kmerlen.try_into().unwrap()); // xor 0xaaaaaaaa and reverse per 2bit
        assert_eq!(0x5F2A_9405_AB86_7BEA, rc, "{:#X}", rc); // TGAC.. 2301.. (in 2bits)

        assert_eq!(xmer.dna.to_usize(), dna); // GTCAGTCAGTCAGTCA => 3210321032103210 (in 2bits)
        assert_eq!(xmer.rc.to_usize(), rc); // xor 0xaaaaaaaa and reverse per 2bit

        let mut mark = xmer.get_hash_and_p(&kc, 0x55);

        // lowest is dna, so that will be hashed by above function.
        let mut test = kc.xmer_hash(dna, 0x55);
        assert_eq!(test, 0xFBFB_3A4A_FABC_021F, "{:#X}", test);

        // its highest 'overbit' is set, so bitwise complementing occurs.
        assert_ne!(test & 0x8000_0000_0000_0000, 0, "{:#x}", test);
        test = !test & 0x7FFF_FFFFFFFF_FFFF;

        assert_eq!(mark.idx, test, "test: {:#x} {:#x}", mark.idx, test);

        let top = (dna & !0x55 & 0xFFFF_FFFF) << 32;
        let bottom = (dna >> 32) & 0x55;
        assert_eq!(
            mark.idx,
            (!dna ^ (top | bottom)) & 0x7FFF_FFFFFFFF_FFFF,
            "top|bottom: {:#x} {:#x}",
            mark.idx,
            (!dna ^ (top | bottom)) & 0x7FFF_FFFFFFFF_FFFF
        );
        assert_eq!(mark.idx, 0x404c5b50543fde0, "{:#x}", mark.idx);
        assert_eq!(mark.p.x(), 0x55);

        let mut undo = mark.idx;
        undo ^= 0xFFFF_FFFFFFFF_FFFF;
        undo = kc.xmer_hash(undo, 0x55);
        assert_eq!(undo, dna, "undo: {:#x} {:#x}", undo, dna);

        kc.extend_xmer(&mut mark).unwrap();
        assert_eq!(mark.p.x(), 0x56);

        let same = xmer.get_hash_and_p(&kc, 0x56);
        assert_eq!(mark.idx, same.idx, "undo: {:#x} {:#x}", mark.idx, same.idx);

        for _ in 0..100 {
            kc.extend_xmer(&mut mark).unwrap();
            let same = xmer.get_hash_and_p(&kc, mark.p.x());
            assert_eq!(mark.idx, same.idx, "undo: {:#x} {:#x}", mark.idx, same.idx);
        }
    }
}
