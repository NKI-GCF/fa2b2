use crate::kmer::xmer_hash;
use crate::kmerloc::ExtPosEtc;
use crate::new_types::position::{BasePos, Position};
use crate::rdbg::STAT_DB;
use num::{FromPrimitive, PrimInt};
use std::cmp;
use std::mem::size_of;

pub struct KmerConst {
    pub no_kmers: usize,
    pub kmerlen: usize,
    pub bitlen: usize,
    pub venster: usize,
    max_afstand: usize,
    pub extent: Vec<usize>,
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
    pub fn from_bitlen(bitlen: usize) -> Self {
        let kmerlen = bitlen / 2;
        let max_afstand = kmerlen / 2;

        // Een venster, groot genoeg voor de meeste unieke x-mers (wat arbitrair)
        let venster = kmerlen + max_afstand;

        // e.g. with a RL 4 & KL 2: (0,1), (1,2), (2,3) => 3 kmers.
        let no_kmers = venster - kmerlen + 1;

        // generate all combinations of 2 kmer positions for extension, high and low bits.
        // overlap and combination duplicates are allowed, those will result respectively
        // in no xor-hash and a complemented xor-hash with the max for index.
        let extent: Vec<usize> = (0..kmerlen).collect();

        dbg_restart!(
            "venster: {}, kmerlen: {}, max_afstand: {}\n--",
            venster,
            kmerlen,
            max_afstand
        );
        if !cfg!(debug_assertions) {
            eprintln!(
                "venster: {}, kmerlen: {}, max_afstand: {}",
                venster, kmerlen, max_afstand
            );
        }
        KmerConst {
            no_kmers,
            kmerlen,
            bitlen,
            venster,
            max_afstand,
            extent,
        }
    }
    pub fn get_ext_max(&self) -> usize {
        self.extent.len()
    }

    // same hash function, from Kmer.
    pub fn get_next_xmer(&self, orig_hash: usize, mut p: ExtPosEtc) -> Option<(usize, ExtPosEtc)> {
        if p.x() + 1 < 0x1_0000_0000 {
            /*was < self.get_ext_max()*/
            let k = self.kmerlen as u32;
            let overbit = 1 << (k * 2 - 1);

            // in get_hash_and_p() bits are flipped if the highest bit was set.
            //
            let mut hash = xmer_hash(orig_hash, p.x(), k);

            if hash < hash.revcmp(k) {
                // XXX: why is this not the inverse ??

                //then flipped, yes: orig_hash here !!
                hash = xmer_hash(!orig_hash & (overbit | (overbit - 1)), p.x(), k);
            }

            p.extend();
            // set to idx for next extension; x is incremented:
            hash = xmer_hash(hash, p.x(), k);
            if (hash & overbit) == 0 {
                Some((hash, p))
            } else {
                Some(((overbit - 1) & !hash, p))
            }
        } else {
            None
        }
    }

    pub fn new(genomesize: usize) -> Self {
        // bit width, required to store all (cumulative) genomic positions, is used as len

        let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
        if (bitlen & 1) == 1 {
            // must be even.
            bitlen += 1
        }
        KmerConst::from_bitlen(bitlen)
    }

    pub fn no_xmers(&self, x: usize) -> usize {
        self.no_kmers - self.afstand(x)
    }
    pub fn get_kmers(&self, x: usize) -> (usize, usize) {
        dbg_assert!(x < self.extent.len());
        let first = self.extent[x] >> self.max_afstand;
        let second = self.extent[x] & (self.max_afstand - 1);
        (first, second)
    }
    pub fn afstand(&self, x: usize) -> usize {
        let kmer = self.get_kmers(x);
        cmp::max(kmer.0, kmer.1)
    }

    pub fn get_kmer_boundaries(
        &self,
        pos: Position,
        contig: (Position, Position),
    ) -> (Position, Position) {
        //let afs = self.afstand(p.x()); // => zou van venster / no_kmers afgetrokken kunnen
        let venster = Position::from(BasePos::from(self.venster));
        let lower = if pos > venster {
            pos - venster
        } else {
            Position::zero()
        };
        (
            cmp::max(lower, contig.0),
            cmp::min(pos + Position::from(BasePos::from(self.no_kmers)), contig.1),
        )
    }
}
