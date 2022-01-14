use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use anyhow::Result;
use num::{FromPrimitive, PrimInt};
use std::cmp;
use std::mem::size_of;

pub struct KmerConst {
    pub(crate) no_kmers: usize,
    pub(crate) kmerlen: usize,
    pub(crate) bitlen: usize,
    pub(crate) read_len: usize,
    pub(crate) dna_topb2_shift: u32,
    pub(crate) rc_mask: u64,
    // put these in kc:
    pub(crate) half_mask: usize,
    pub(crate) overbit: usize,
    pub(crate) xmer_mask: usize,
    pub(crate) over_mask: usize,
    // XXX!seed influences order of selection, a good seed may improve indexing / mapping.
    // XXX!HOWEVER!XXX indexing and mapping should use the same seed !!
    pub(crate) seed: usize,
}

pub trait RevCmp<T: PrimInt + FromPrimitive> {
    fn revcmp(self, kmerlen: u32) -> T;
}

/// create bitmask. e.g. dvm::<u32>(0xf0  0xff) => 0xf0_f0_f0_f0
#[inline]
fn dvm<T: PrimInt + FromPrimitive>(numerator: u32, divisor: u32) -> T {
    let base = T::max_value() / T::from_u32(divisor).unwrap();
    T::from_u32(numerator).unwrap() * base
}

#[macro_export]
macro_rules! implement_revcmp { ($($ty:ty),*) => ($(
    /// give twobit reverse complent for given kmerlen
    impl RevCmp<$ty> for $ty {
        #[inline]
        fn revcmp(self, kmerlen: u32) -> $ty {
            let mut seq = self.swap_bytes() ^ dvm::<$ty>(2, 3);

            let left_nibbles = (seq & dvm::<$ty>(0xf, 0xff)) << 4;
            let right_nibbles = (seq & dvm::<$ty>(0xf0, 0xff)) >> 4;
            seq = left_nibbles | right_nibbles;

            let all_left_two_bits = (seq & dvm::<$ty>(0x3, 0xf)) << 2;
            let all_right_two_bits = (seq & dvm::<$ty>(0xc, 0xf)) >> 2;
            seq = all_left_two_bits | all_right_two_bits;

            seq >> (size_of::<$ty>() * 8 - usize::try_from(kmerlen).unwrap() * 2)
        }
    }
    )*)
}
implement_revcmp!(usize, u64);

impl KmerConst {
    pub(crate) fn new(genomesize: usize, read_len: usize) -> Self {
        // bit width, required to store all (cumulative) genomic positions, is used as len

        let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
        if (bitlen & 1) == 1 {
            // must be even.
            bitlen += 1
        }
        KmerConst::from_bitlen(bitlen, read_len)
    }

    pub(crate) fn from_bitlen(bitlen: usize, read_len: usize) -> Self {
        let kmerlen = bitlen / 2;
        dbg_assert!(kmerlen > 0);
        dbg_assert!(read_len >= kmerlen);

        // e.g. with a RL 4 & KL 2: (0,1), (1,2), (2,3) => 3 kmers.
        let no_kmers = read_len + 1 - kmerlen;

        dbg_restart!("read_len: {}, kmerlen: {}\n--", read_len, kmerlen);
        if !cfg!(debug_assertions) {
            eprintln!("read_len: {}, kmerlen: {}", read_len, kmerlen);
        }
        let dna_topb2_shift = (u32::try_from(kmerlen).expect("kmer too large") << 1) - 2;
        let overbit = 1_usize
            .checked_shl(dna_topb2_shift + 1)
            .expect("k-mer shift");
        let xmer_mask = overbit - 1;
        KmerConst {
            no_kmers,
            kmerlen,
            bitlen,
            read_len,
            dna_topb2_shift,
            rc_mask: (1 << dna_topb2_shift) - 1,
            half_mask: (1 << kmerlen) - 1,
            overbit,
            xmer_mask,
            over_mask: overbit | xmer_mask,
            seed: 0x5EED,
        }
    }
    #[inline(always)]
    pub(crate) fn xmer_hash(&self, idx: usize, seed: usize) -> usize {
        idx ^ (((idx & !seed & self.half_mask) << self.kmerlen) | ((idx >> self.kmerlen) & seed))
    }
    pub(crate) fn hash_and_compress(&self, seq: usize, x: usize) -> usize {
        self.compress_xmer(self.xmer_hash(seq, x))
    }

    // same hash function, from Xmer.
    pub(crate) fn get_next_xmer(
        &self,
        orig_hash: usize,
        mut p: ExtPosEtc,
    ) -> Result<(usize, ExtPosEtc)> {
        let old_x = p.x();
        p.extend()?;

        // in get_hash_and_p() bits are flipped if the highest bit was set.
        let mut hash = self.xmer_hash(orig_hash, old_x);

        if hash < hash.revcmp(self.kmerlen as u32) {
            // XXX: why is this not the inverse ??

            //then flipped, yes: orig_hash here !!
            hash = self.xmer_hash(self.over_mask & !orig_hash, old_x);
        }

        // set to idx for next extension; x is incremented:
        hash = self.xmer_hash(hash, p.x());
        if (hash & self.overbit) == 0 {
            Ok((hash, p))
        } else {
            Ok((self.xmer_mask & !hash, p))
        }
    }
    /// compress base_seq or xmer_hash from xmer module
    pub(crate) fn compress_xmer(&self, v: usize) -> usize {
        dbg_assert!(v & self.over_mask == v);
        if v & self.overbit == 0 {
            v
        } else {
            self.xmer_mask & !v
        }
    }

    pub(crate) fn get_kmer_boundaries(
        &self,
        pos: Position,
        contig: (Position, Position),
    ) -> (Position, Position) {
        let read_len = Position::from(BasePos::from(self.read_len));
        let lower = if pos > read_len {
            pos - read_len
        } else {
            Position::zero()
        };
        (
            cmp::max(lower, contig.0),
            cmp::min(pos + Position::from(BasePos::from(self.no_kmers)), contig.1),
        )
    }
}
