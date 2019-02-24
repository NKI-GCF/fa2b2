use std::mem::size_of;
use bit_reverse::ParallelReverse;
use kmerloc::PriExtPosOri;

pub struct KmerConst {
	pub max_no_kmers: usize,
	pub kmerlen: usize,
	pub bitlen: usize,
	pub readlen: usize,
	pub ext_max: usize,
}

pub fn afstand(x: usize, kmerlen: usize) -> usize {
	let t = 1 << x;
	if t <= kmerlen {t} else {
		let n = kmerlen.next_power_of_two().trailing_zeros() as usize;
		n + (x - n) * kmerlen
	}
}

impl KmerConst {
	pub fn new(readlen: usize, genomesize: usize) -> Self {
		// bit width, required to store all (cumulative) genomic positions, is used as len
		println!("Genome size: {}", genomesize);

		let mut bitlen = genomesize.next_power_of_two().trailing_zeros() as usize;
		if (bitlen & 1) == 1 { // must be even.
			bitlen += 1
		}
		let kmerlen = bitlen / 2;
		let max_no_kmers = readlen - kmerlen;
		let mut p = 0;
		while {
			p.extend();
			afstand(p.x(), kmerlen) <= max_no_kmers
		}{}
		println!("Using a kmerlength of {}, readlength of {}, ext_max: {}\n--", bitlen / 2, readlen, p.x());
		KmerConst {
			max_no_kmers,
			kmerlen,
			bitlen,
			readlen,
			ext_max: p.x()
		}
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
