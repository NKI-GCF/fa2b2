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
    pub fn new(genomesize: usize) -> Self {
        // bit width, required to store all (cumulative) genomic positions, is used as len

        let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
        if (bitlen & 1) == 1 {
            // must be even.
            bitlen += 1
        }
        let kmerlen = bitlen / 2;
        let max_afstand = kmerlen / 2;
        let venster = kmerlen + max_afstand;

        // e.g. with a RL 4 & KL 2: (0,1), (1,2), (2,3) => 3 kmers.
        let no_kmers = venster - kmerlen + 1;

        // generate all combinations of 2 kmer positions for extension, high and low bits.
        // overlap and combination duplicates are allowed, those will result respectively
        // in no xor-hash and a complemented xor-hash with the max for index.
        let extent: Vec<usize> = (0..kmerlen).collect();

        dbg_restart!(
            "Genome size: {}, venster: {}, kmerlen: {}, max_afstand: {}\n--",
            genomesize,
            venster,
            kmerlen,
            max_afstand
        );
        if !cfg!(debug_assertions) {
            eprintln!(
                "Genome size: {}, venster: {}, kmerlen: {}, max_afstand: {}",
                genomesize, venster, kmerlen, max_afstand
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
    pub fn get_kmers(&self, x: usize) -> (usize, usize) {
        let first = self.extent[x] >> self.max_afstand;
        let second = self.extent[x] & (self.max_afstand - 1);
        (first, second)
    }
    pub fn afstand(&self, x: usize) -> usize {
        let kmer = self.get_kmers(x);
        cmp::max(kmer.0, kmer.1)
    }

    pub fn get_kmer_boundaries(&self, p: u64, contig: (u64, u64)) -> (u64, u64) {
        (
            cmp::max(p.saturating_sub(self.venster as u64 * 2), contig.0),
            cmp::min(p + (self.no_kmers as u64 - 1) * 2, contig.1),
        )
    }
}
