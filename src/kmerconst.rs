use crate::kmerloc::PriExtPosOri;
use crate::rdbg::STAT_DB;
use std::cmp;

pub struct KmerConst {
    pub no_kmers: usize,
    pub kmerlen: usize,
    pub bitlen: usize,
    pub readlen: usize,
    pub ext_max: usize,
    pub extension: [(u8, u8); 256],
}

impl KmerConst {
    pub fn new(readlen: usize, genomesize: usize) -> Self {
        // bit width, required to store all (cumulative) genomic positions, is used as len

        let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
        if (bitlen & 1) == 1 {
            // must be even.
            bitlen += 1
        }
        let kmerlen = bitlen / 2;
        // e.g. with a RL 4 & KL 2: (0,1), (1,2), (2,3) => 3 kmers.
        let no_kmers = readlen - kmerlen + 1;

        // generate all combinations of 2 kmer positions for extension. overlap and inverse
        // combinations are allowed, those will result respectively in no xor-hash and a
        // complemented xor-hash with the max for index.
        // the order used here maintains the minimum (positional) distance betweens kmers
        // for all extensions
        let mut extension = [(0, 0); 256];
        let mut i = 0;
        let mut j = 0;
        for element in extension.as_mut_slice() {
            *element = (i, j);
            match i.cmp(&j) {
                cmp::Ordering::Less => i += 1,
                cmp::Ordering::Greater => {
                    j += 1;
                    if i == j {
                        i = 0;
                    }
                }
                cmp::Ordering::Equal => {
                    j = 0;
                    i += 1;
                }
            }
        }

        dbg_restart!(
            "Genome size: {}, readlen: {}, kmerlen: {}, ext_max: {}\n--",
            genomesize,
            readlen,
            kmerlen,
            0xff
        );
        KmerConst {
            no_kmers,
            kmerlen,
            bitlen,
            readlen,
            ext_max: 0xff,
            extension,
        }
    }
    pub fn afstand(&self, x: usize) -> usize {
        cmp::min(
            cmp::max(self.extension[x].0, self.extension[x].1) as usize,
            self.readlen,
        )
    }

    pub fn get_kmer_boundaries(&self, p: u64, contig: (u64, u64)) -> (u64, u64) {
        (
            cmp::max(p.saturating_sub(self.readlen as u64 * 2), contig.0),
            cmp::min(p + (self.no_kmers as u64 - 1) * 2, contig.1),
        )
    }

    pub fn leftmost_of_scope(&self, p: u64, plim_0: u64) -> u64 {
        cmp::max(
            plim_0 | p.extension(),
            p.saturating_sub((self.no_kmers + self.afstand(p.x())) as u64),
        )
    }
}
