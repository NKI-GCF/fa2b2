use crate::kmerloc::{BasePos, ExtPosEtc};
use crate::new_types::position::Position;
use crate::rdbg::STAT_DB;
use std::cmp;

pub struct KmerConst {
    pub no_kmers: usize,
    pub kmerlen: usize,
    pub bitlen: usize,
    pub venster: usize,
    max_afstand: usize,
    pub extent: Vec<usize>,
}

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
            cmp::min(pos + (self.no_kmers as u64 - 1).basepos_to_pos(), contig.1),
        )
    }
}
