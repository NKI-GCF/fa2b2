use crate::kmerloc::PriExtPosOri;
use crate::rdbg::STAT_DB;
use bit_reverse::ParallelReverse;
use std::cmp::min;
use std::mem::size_of;

pub struct KmerConst {
    pub no_kmers: usize,
    pub kmerlen: usize,
    pub bitlen: usize,
    pub readlen: usize,
    pub ext_max: usize,
}

fn afstand(x: usize, kmerlen: usize) -> usize {
    let t = if x == 0 { 0 } else { 1 << (x - 1) };
    if t <= kmerlen {
        t
    } else {
        let n = kmerlen.next_power_of_two().trailing_zeros() as usize;
        // (x - log2(kmerlen)) * kmerlen
        (x - n) * kmerlen
    }
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
        let mut p = 0;
        while {
            p.extend();
            afstand(p.x(), kmerlen) < no_kmers
        } {}
        dbg_restart!(
            "Genome size: {}, readlen: {}, kmerlen: {}, ext_max: {}\n--",
            genomesize,
            readlen,
            kmerlen,
            p.x()
        );
        KmerConst {
            no_kmers,
            kmerlen,
            bitlen,
            readlen,
            ext_max: p.x(),
        }
    }
    pub fn afstand(&self, x: usize) -> usize {
        min(afstand(x, self.kmerlen), self.readlen)
    }

    pub fn get_kmer_boundaries(&self, p: u64, contig_start: u64, contig_end: u64) -> (u64, u64) {
        let p_rl = (self.readlen << 1) as u64;

        let left = if p >= contig_start + p_rl {
            p - p_rl
        } else {
            contig_start
        };
        let right = min(p + ((self.readlen - self.kmerlen) << 1) as u64, contig_end);
        (left, right)
    }

    /// with given extension, create bitmask to flip high bits before extreme minimization
    /// with this each extension minimizes in its own domain, decreasing ks.kmp sparsity.
    // one added because index is shortened (kmer index top bit flipped, if set)
    // XXX could use a lookup array instead?
    pub fn ext_domain(&self, x: usize) -> usize {
        x.swap_bits() >> ((size_of::<usize>() << 3) - self.bitlen + 1)
    }
}
