use std::mem::size_of;
use bit_reverse::ParallelReverse;

pub struct KmerConst {
	pub pos_mask: u64,
	pub max_no_kmers: usize,
	pub kmerlen: usize,
	pub bitlen: usize,
	pub readlen: usize,
	pub priority_shft: usize,
	pub extent: u32,
}

impl KmerConst {
	pub fn new(readlen: usize, genomesize: usize, extent: u32) -> Self {
		// bit width, required to store all (cumulative) genomic positions, is used as len
		println!("Genome size: {}", genomesize);

		let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
		if (bitlen & 1) == 1 {
			bitlen += 1
		}
		let kmerlen = bitlen / 2;
		let max_no_kmers = readlen - kmerlen;
		println!("Using a kmerlength of {}, readlength of {}, ext_max: {}\n--", bitlen / 2, readlen, (max_no_kmers + 1).next_power_of_two().trailing_zeros());
		KmerConst {
			pos_mask: (1 << extent) - 1,
			max_no_kmers,
			kmerlen,
			bitlen,
			readlen,
			priority_shft: 8 * size_of::<u64>() - 1,
			extent,
		}
	}
	pub fn priority(&self) -> u64 {
		1 << self.priority_shft
	}
	pub fn ext_max(&self) -> usize {
		//(size_of::<usize>() * 8) - (self.max_no_kmers.leading_zeros() as usize) - 1
		(self.max_no_kmers + 1).next_power_of_two().trailing_zeros() as usize
	}

	/// with given extension, create bitmask to flip high bits before extreme minimization
	/// with this each extension minimizes in its own domain, decreasing ks.kmp sparsity.
	// one added because index is shortened (kmer index top bit flipped, if set)
	pub fn ext_domain(&self, x: usize) -> usize {
		x.swap_bits() >> (size_of::<usize>() * 8 - self.bitlen + 1)
	}
}

#[cfg(test)]
mod tests {
    const READLEN: usize = 64;
    const KMERLEN: usize = 16;
    #[test]
    fn ext_max() {
        assert_eq!(1, 1);
    }
}
