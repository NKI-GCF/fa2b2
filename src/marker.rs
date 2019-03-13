extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate num_traits;
extern crate num;

use std::cmp;
use std::slice::Iter;

use occurrence::Occurrence;
use kmerloc::{PriExtPosOri,MidPos};
use kmerstore::KmerStore;

pub struct KmerIter<'a> {
	n_stretch: u64,
	goffs: u64,
	pub(super) occ: &'a mut Vec<Occurrence<'a>>,
	pub(super) ks: &'a mut KmerStore<u64>,
} //^-^\\

impl<'a> KmerIter<'a> {
	pub fn new(ks: &'a mut KmerStore<u64>, occ: &'a mut Vec<Occurrence<'a>>) -> Self {
		KmerIter { n_stretch: 0, goffs: 0, occ, ks }
	}

	fn finalize_n_stretch_if_pending(&mut self) {
		if self.n_stretch > 0 {
			eprintln!("added new contig. Ns:{}", self.n_stretch);
			self.ks.offset_contig(self.n_stretch);
			self.goffs += self.n_stretch;
			self.n_stretch = 0;
		}
	}

	/// if N, insert contig (once) and count stretch, false
	/// else store and update occurrance kmers, complete ? true : false
	fn complete_occurance_or_contig(
		&mut self,
		n: usize,
		b2: u8,
	) -> bool {
		if let Some(occ) = self.occ.get_mut(n) {
			let p = occ.p.pos();
			if n != 0 {
				debug_assert!((self.ks.b2[p as usize >> 3] >> (p & 6)) & 3 == b2);
				return occ.complete(self.ks, b2, n);
			}
			if b2 < 4 {
				if let Some(qb) = self.ks.b2.get_mut(p as usize >> 3) {
					*qb |= b2 << (p & 6);
				}
				self.ks.p_max = (p + 2) & !1;
				if occ.complete(self.ks, b2, n) {
					return true;
				}
			} else if occ.i != 0 {
				self.goffs += occ.i as u64;
				eprintln!("started N-stretch at {}.", p);
				self.ks.push_contig(p, self.goffs);

				// clear all except orientation and position to rebuild at the start of a new contig.
				//assert_eq!(occ.i, 1);
				occ.i = 0;
				self.n_stretch = 1;
				return false;
			} else {
				self.n_stretch += 1;
				return false;
			}
		}
		self.finalize_n_stretch_if_pending();
		false
	}

	/// when rebuilding eq for recurrent kmer, and extending take into account contig boundaries
	/// for that site
	fn get_plimits(&self, p: u64) -> (u64, u64) {

		// binary search; limit endp to end of contig
		debug_assert!(p != 0);
		let i = self.ks.get_contig(p);

		debug_assert!(i < self.ks.contig.len());
		debug_assert!(self.ks.contig[i].twobit <= p);

		let contig_start = self.ks.get_twobit_before(i).unwrap_or(0);
		let kc = self.occ[0].kc;
		let p_max = self.ks.p_max;

		let p_no_kmers = ((kc.no_kmers - kc.kmerlen) << 1) as u64;

		let left = if dbgx!(p >= contig_start + p_no_kmers) {
			p - p_no_kmers
		} else {
			contig_start
		};

		let right = cmp::min(
			p + (kc.readlen << 1) as u64,
			self.ks.get_twobit_after(i).unwrap_or(p_max),
		);
		(left, right)
	}

	fn rebuild_kmer_stack(&self, stored_at_index: &u64) -> Occurrence<'a> {

		let plimits = self.get_plimits(stored_at_index.pos());
		Occurrence::new(plimits, self.occ[0].kc, stored_at_index.extension())
	}

	fn search_occ_for_pos(&self, original_n: usize, stored_at_index: u64) -> usize {


		for n in 0..original_n {
			if dbgf!(stored_at_index == self.occ[n].mark.p,
					"{:#?}: {:#016x} == {:#016x}?", stored_at_index, self.occ[n].mark.p) {
				return dbgx!(n);
			}
		}
		original_n
	}
	fn add_newstack(&mut self, next_stack: Occurrence<'a>, n: usize) {
		if n < self.occ.len() {
			self.occ[n] = next_stack;
		} else {
			self.occ.push(next_stack);
		}
	}

	/// when the stored entry in ks.kmp is replaced or blacklisted, the occurence for the stored
	/// needs to be extended to find, there, a minposfor the extended kmers.
	/// the entry is added to the stack, return the index.
	fn get_next_for_extension(&mut self, stored_at_index: &mut u64, original_n: usize,
						   is_replaced: bool) -> usize {

		let mut n = self.search_occ_for_pos(original_n, *stored_at_index);
		if dbgx!(*stored_at_index != self.occ[n].mark.p) {

			stored_at_index.extend();
			let next_stack = self.rebuild_kmer_stack(dbgf!(stored_at_index, "{:#x}"));

			if is_replaced && n != 0 {
				// position is written below, this self.occ[n] is done.
				// Could add it to recurring kmer_stacks, though.
				self.occ[n] = next_stack;
			} else {
				// 0th is never overwritten.
				n += 1;
				self.add_newstack(next_stack, n);
			}
		}
		n
	}

	fn next_b2(&self, seq: &mut Iter<u8>, n: usize) -> Option<u8> {
		if n == 0 {
			seq.next().map(|c| (c >> 1) & 7u8)
		} else {
			self.ks.b2_for_p(self.occ[n].p.pos())
		}
	}

	pub fn markcontig<T: MidPos>(
		&mut self,
		seq: &mut Iter<u8>,
	) -> u64 {

		let mut n = 0;
		self.goffs = 0;
		self.ks.push_contig(self.occ[n].p.pos(), self.goffs);

		while let Some(b2) = self.next_b2(seq, n) {

			if !self.complete_occurance_or_contig(n, b2) {
				continue;
			}
			loop {

				let (min_index, min_pos) = (self.occ[n].mark.idx, self.occ[n].mark.p);

				let mut stored_at_index = self.ks.kmp[min_index];

				if dbgx!(stored_at_index.is_replaceable_by(min_pos)) {
					eprintln!("[{:#x}] (={:#016x}) <= {:#016x}", min_index, stored_at_index, min_pos);

					self.ks.kmp[min_index] = min_pos;

					if dbgx!(stored_at_index.is_set_and_not(min_pos)) {

						n = self.get_next_for_extension(&mut stored_at_index, n, true);

					} else if n > 0 {
						n -= 1;
						continue;
					}

					break; // position written (done) or added next_stack which requires extension.
				}
				if dbgx!(self.occ[n].try_extension_redefine_minimum()) {

					if dbgx!(stored_at_index.extension() == min_pos.extension()) {
						self.ks.kmp[min_index].blacklist();

						// both stored and current require extensions
						n = self.get_next_for_extension(&mut stored_at_index, n, false);

						break; // next_stack .
					}
				} else {
					self.ks.kmp[min_index].blacklist();
					if dbgx!(n == 0) { // never pop 0th. 0th needs to be renewed when done.
						break;
					}
					n -= 1;
					if dbgx!(n == 0 && self.occ[n].mark.p == self.ks.kmp[self.occ[n].mark.idx]) {
						// test already done
						break;
					}
				}
			}
		}
		self.finalize_n_stretch_if_pending();
		self.occ[0].p.pos() as u64
	}
}

#[cfg(test)]
mod tests {
	use kmerconst::KmerConst;
	use super::KmerIter;
	use super::KmerStore;
	use super::Occurrence;
	use super::PriExtPosOri;
	const READLEN: usize = 16;
	const SEQLEN: usize = 250;

	fn process<'a>(occ: &'a mut Vec<Occurrence<'a>>, ks: &'a mut KmerStore<u64>, seq: Vec<u8>) {
		let mut kmi = KmerIter::new(ks, occ);
		kmi.markcontig::<u64>(&mut seq.iter());
	}

	#[test]
	fn test_16n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NNNNNNNNNNNNNNNN"[..].to_owned());
		}
		assert_eq!(ks.contig.len(), 1);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 16);
	}
	#[test]
	fn test_1n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"N"[..].to_owned());
		}
		assert_eq!(ks.contig.len(), 1);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 1);
	}
	#[test]
	fn test_1n1c1n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NCN"[..].to_owned());
		}
		assert_eq!(ks.contig.len(), 2);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 1);
		assert_eq!(ks.contig[1].twobit, 2);
		assert_eq!(ks.contig[1].genomic, 3);

	}
	#[test]
	fn test_17c() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"CCCCCCCCCCCCCCCCC"[..].to_owned());
		}
		assert_eq!(ks.kmp.len(), 128);
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let unresolveable = (1 << 63) | ((kc.ext_max as u64) << 48);
		assert_eq!(ks.kmp[0], unresolveable, "test_17_c(): {:x}", ks.kmp[0]);
	}
	#[test]
	fn test_1n18c1n() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NCCCCCCCCCCCCCCCCCCN"[..].to_owned());
		}
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let unresolveable = (1 << 63) | ((kc.ext_max as u64) << 48);
		assert_eq!(ks.kmp[0], unresolveable);
	}
	#[test]
	fn test_1n16c() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"NCCCCCCCCCCCCCCCC"[..].to_owned());
		}
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let first_pos = (1 << 63) | ((kc.kmerlen as u64)<< 1);
		assert_eq!(ks.kmp[0], first_pos);
	}
	#[test]
	fn test_16at1a() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			process(&mut occ, &mut ks, b"ATATATATATATATATA"[..].to_owned());
		}
		for i in 0..ks.kmp.len() {
			if i != 0x0 && i != 0x22 {
				assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
			}
		}
		let unresolveable = (1 << 63) | ((kc.ext_max as u64) << 48);
		assert_eq!(ks.kmp[0], unresolveable);
	}
	#[test]
	fn test_reconstruct1() {
		let kc = KmerConst::new(READLEN, SEQLEN);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		let ks_kmp_len = ks.kmp.len();
		let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
		let mut kmi = KmerIter::new(&mut ks, &mut occ);
		let seq_vec =
			b"GCGATATTCTAACCACGATATGCGTACAGTTATATTACAGACATTCGTGTGCAATAGAGATATCTACCCC"[..].to_owned();
		let mut seq = seq_vec.iter();
		kmi.markcontig::<u64>(&mut seq);
		for hash in 0..ks_kmp_len {
			let mut p = kmi.ks.kmp[hash];
			if !p.blacklisted() {
				let new_stack = kmi.rebuild_kmer_stack(&p);
				kmi.add_newstack(new_stack, 1);
				while let Some(b2) = kmi.next_b2(&mut seq, 1) {
					if kmi.complete_occurance_or_contig(1, b2) {
						break;
					}
				}
				let mark_p = kmi.occ[1].mark.p;
				eprintln!("testing: [{:#x}]: {:#x} == (stored p:){:#x}", hash, mark_p, p);
				assert_eq!(mark_p, p);
			}
		}
	}
	#[test]
	fn test_reconstruct_simple() { // all mappable.
		let kc = KmerConst::new(4, 15);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		let ks_kmp_len = ks.kmp.len();
		let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
		let mut kmi = KmerIter::new(&mut ks, &mut occ);
		let seq_vec = b"GGAACCTTCAGAGTG"[..].to_owned();
		let mut seq = seq_vec.iter();
		kmi.markcontig::<u64>(&mut seq);
		for hash in 0..ks_kmp_len {
			let mut p = kmi.ks.kmp[hash];
			if !p.blacklisted() {
				let new_stack = kmi.rebuild_kmer_stack(&p);
				kmi.add_newstack(new_stack, 1);
				while let Some(b2) = kmi.next_b2(&mut seq, 1) {
					if kmi.complete_occurance_or_contig(1, b2) {
						break;
					}
				}
				let mark_p = kmi.occ[1].mark.p;
				eprintln!("testing: [{:#x}]: {:#x} == (stored p:){:#x}", hash, mark_p, p);
				assert_eq!(mark_p, p);
			}
		}
	}
	#[test]
	fn test_reconstruct_simplest() { // all mappable.

		let mut errors = 0;
		for gen in 0..=255 {
			let kc = KmerConst::new(2, 4);
			let mut ks = KmerStore::<u64>::new(kc.bitlen);
			let ks_kmp_len = ks.kmp.len();
			let mut occ: Vec<Occurrence> = vec![Occurrence::new((0, u64::max_value()), &kc, 0)];
			let mut kmi = KmerIter::new(&mut ks, &mut occ);
			let seq_vec:Vec<_> = [gen & 3, (gen >> 2) & 3, (gen >> 4) & 3, (gen >> 6) & 3].iter().map(|b| match b {
					0 => 'A', 1 => 'C', 2 => 'T', 3 => 'G', _ => panic!("here") }).collect();
			eprint!("\nsequence: {:?}", seq_vec);

			eprint!("\n");
			let vv: Vec<u8> = seq_vec.iter().map(|&c| c as u8).collect();
			let mut seq = vv.iter();
			kmi.markcontig::<u64>(&mut seq);
			let mut miss = 0;
			for hash in 0..ks_kmp_len {
				let mut p = kmi.ks.kmp[hash];
				if !p.blacklisted() {
					//dbgf!(hash, "[{:#x}] = {:#x}", p);
					let new_stack = kmi.rebuild_kmer_stack(&p);
					kmi.add_newstack(new_stack, 1);
					eprintln!("occ.p={:x} .plim={:x}", kmi.occ[1].p, kmi.occ[1].plim);
					while let Some(b2) = kmi.next_b2(&mut seq, 1) {
						eprintln!("occ.p={:x} twobit {:x}", kmi.occ[1].p, b2);
						if kmi.complete_occurance_or_contig(1, b2) {
							break;
						}
					}
					let mark_p = kmi.occ[1].mark.p;
					eprintln!("testing: [{:#x}]: {:#x} == (stored p:){:#x}", hash, mark_p, p);
					if mark_p != p {
						miss += 1;
					}
				}
			}
			if miss != 0 {
				eprintln!("observed {} incorrect rebuilds", miss);
				errors += 1;
			}
		}
		eprintln!("observed errors in {}/256", errors);
		assert!(false);
	}
}
