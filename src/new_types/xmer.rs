use super::{
    extension::Extension,
    position::Position,
    twobit::{TwoBit, TwoBitDna, TwoBitRcDna},
};
use crate::new_types::extended_position::ExtPosEtc;
use crate::rdbg::STAT_DB;
use derive_more::Into;
use std::cmp;

#[derive(Copy, Clone, PartialEq, Eq)]
/// A kmer that dissociates index and strand orientation
pub struct Xmer {
    pub(super) dna: TwoBitDna,
    pub(super) rc: TwoBitRcDna,
    pub pos: Position,
    half_mask: usize,
    overbit: usize,
    xmer_mask: usize,
    pub kmerlen: u32,
} //^-^\\

#[inline(always)]
pub(crate) fn xmer_hash(idx: usize, x: usize, half_mask: usize, k: u32) -> usize {
    idx ^ (((idx & !x & half_mask) << k) | ((idx >> k) & x))
}

impl Xmer {
    /// get a kmer for this length
    pub(crate) fn new(kmerlen: u32) -> Self {
        let overbit = 1_usize << (kmerlen * 2 - 1);
        Xmer {
            dna: TwoBitDna::new(0),
            rc: TwoBitRcDna::new(0),
            pos: Position::zero(),
            half_mask: (1 << kmerlen) - 1,
            overbit,
            xmer_mask: overbit - 1,
            kmerlen,
        }
    }

    /// true if the kmer is from the template. Palindromes are special.
    #[inline(always)]
    pub(crate) fn is_template(&self) -> bool {
        self.dna.lt_strand(self.rc)
    }

    /// adds twobit to k-mer sequences, to dna in the top two bits. Returns orientation.
    pub(crate) fn update(&mut self, b2: TwoBit) -> bool {
        // XXX function is hot
        self.dna.add(b2, self.kmerlen * 2 - 2);
        self.rc.add(b2, self.kmerlen * 2 - 2);
        self.pos.incr();
        self.dna.lt_strand(self.rc)
    }
    //TODO: remove use_min
    #[inline(always)]
    fn get_base_seq(&self) -> usize {
        if self.is_template() {
            self.dna.to_usize()
        } else {
            self.rc.to_usize()
        }
    }
    /// an extension specific index
    pub(crate) fn get_hash(&self, x: usize) -> usize {
        let seq = self.get_base_seq();

        let idx = xmer_hash(seq, x, self.half_mask, self.kmerlen);

        if idx & self.overbit == 0 {
            idx
        } else {
            self.xmer_mask & !idx
        }
    }
    pub(crate) fn get_hash_and_p(&self, x: usize) -> (usize, ExtPosEtc) {
        (
            self.get_hash(x),
            ExtPosEtc::from((Extension::from(x), self.pos)),
        )
    }

    /// a sequence specific index, the same when read in other orientation
    pub(super) fn get_baseidx(&self) -> BaseIdx {
        let seq = self.get_base_seq();

        // flipped if the top bit is set, to reduce size.
        if (seq & self.overbit) == 0 {
            BaseIdx(seq)
        } else {
            BaseIdx(self.xmer_mask & !seq)
        }
    }
}

impl PartialOrd for Xmer {
    fn partial_cmp(&self, other: &Xmer) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Into)]
pub(super) struct BaseIdx(usize);

impl Ord for Xmer {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        self.get_baseidx().0.cmp(&other.get_baseidx().0)
    }
}
