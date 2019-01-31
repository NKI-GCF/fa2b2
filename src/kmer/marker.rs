extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate kmerconst;
extern crate kmerloc;
extern crate kmerstore;
extern crate occurrence;
extern crate num_traits;
extern crate num;

use self::kmerconst::KmerConst;
use self::kmerloc::{PriExtPosOri,MidPos};
use self::kmerstore::KmerStore;
use self::occurrence::Occurrence;
use std::cmp;
use std::slice::Iter;

#[macro_export]
macro_rules! dbgf {
	($l:literal) => ({
		eprint!("[{}:{}] {}\n", file!(), line!(), stringify!($l));
		assert!(cfg!(debug_assertions));
		$l
	});
	($fmt:literal, $expr:expr$(, $opt:expr)*) => {
		match $expr {
			expr => {
				eprint!(concat!("[{}:{}] {} = ", $fmt, "\n"), file!(), line!(), stringify!($expr), &expr$(, $opt)*);
				assert!(cfg!(debug_assertions));
				expr
			}
		}
	}
}



pub struct KmerIter {
	n_stretch: u64,
	goffs: u64,
} //^-^\\

impl<'a> KmerIter {
	pub fn new() -> Self {
		KmerIter { n_stretch: 0, goffs: 0 }
	}

	fn finalize_n_stretch_if_pending(&mut self, ks: &mut KmerStore<u64>) {
		if self.n_stretch > 0 {
			eprintln!("added new contig. Ns:{}", self.n_stretch);
			ks.offset_contig(self.n_stretch);
			self.goffs += self.n_stretch;
			self.n_stretch = 0;
		}
	}

	/// if N, insert contig (once) and count stretch, false
	/// else store and update occurrance kmers, complete ? true : false
	fn complete_occurance_or_contig(
		&mut self, ks: &mut KmerStore<u64>,
		occ: &mut Occurrence,
		b2: u8,
	) -> bool {
		let p = occ.p.pos();
		if b2 < 4 {
			if let Some(qb) = ks.b2.get_mut(p as usize >> 3) {
				*qb |= b2 << (p & 6);
			}
			ks.p_max = (p + 2) & !1;
			if occ.complete(b2) {// ks, 
				return true;
			}
			self.finalize_n_stretch_if_pending(ks);
		} else if occ.i != 0 {
			self.goffs += occ.i as u64;
			eprintln!("started N-stretch at {}.", p);
			ks.push_contig(p, self.goffs);

			// clear all except orientation and position to rebuild at the start of a new contig.
			//assert_eq!(occ.i, 1);
			occ.d.clear();
			occ.i = 0;
			self.n_stretch = 1;
		} else {
			self.n_stretch += 1;
		}
		false
	}

	/// when rebuilding eq for recurrent kmer, and extending take into account contig boundaries
	/// for that site
	fn get_plimits(&self, ks: &KmerStore<u64>, kc: &KmerConst, p: u64) -> (u64, u64) {
		// binary search; limit endp to end of contig
		debug_assert!(p != 0);
		let i = ks.get_contig(p);
		debug_assert!(i < ks.contig.len());
		//println!("contig[{}].twobit:{}, p:{}", i, ks.contig[i].twobit, p);

		debug_assert!(ks.contig[i].twobit <= p);
		let x = ks.get_twobit_before(i).unwrap_or(0);
		let addl = ((kc.max_no_kmers + kc.kmerlen) << 1) as u64;
		let p_max = ks.p_max;
		//eprintln!("p:{:x} > x:{:x} + addl:{:x} ?",  p, x, addl);
		let left = if p > x + addl { p - addl } else { x };
		let right = cmp::min(
			p + (kc.readlen << 1) as u64,
			ks.get_twobit_after(i).unwrap_or(p_max),
		);
		(left, right)
	}

	pub fn markcontig<T: MidPos>(
		&mut self, ks: &'a mut KmerStore<u64>,
		seq: &mut Iter<u8>,
		kc: &'a KmerConst,
	) -> u64 {
		let mut occ = vec![Occurrence::new((0, u64::max_value()), kc, 0, true)];
		let mut q = 0;
		self.goffs = 0;
		ks.push_contig(occ[0].p.pos(), self.goffs);
		while let Some(c) = seq.next() {
			if !self.complete_occurance_or_contig(ks, &mut occ[0], (*c >> 1) & 0x7) {
				continue;
			}
			loop {
				eprintln!("---\n");
				if !occ[q].mark.is_set() {
					let x = occ[q].p.x();
					let offs = if x == 0 {0} else {1<<x};
					let end_i = kc.max_no_kmers - offs;
					for i in 0..end_i {
						let (hash, p) = occ[q].get_hash_p(i, x, (end_i + offs - i) << 1);
						if dbg!(occ[q].hash_is_extreme(hash, p, x)) {
							break;
						}
					}
					dbgf!("{:#x?}", occ[q].mark.p);
				}

				let mark_idx = occ[q].mark.idx;
				let mark_p = occ[q].mark.p;
				let mut former = ks.kmp[mark_idx];
				debug_assert!(mark_p != 0);
				//eprintln!("{:#x?}", former_stack);

				if dbgf!("{:#?}: had stored kmer already this position?", former.is_same(mark_p | (1 << 63))) {
					break;
				} else if dbgf!("{:#?}: former:{:#016x}, mark_p:{:#016x}", former.is_replaceable_by(mark_p), former, mark_p) {

					ks.kmp[mark_idx] = dbgf!("{:x}: position set for kmer [{:#x}]", mark_p, mark_idx);
					ks.kmp[mark_idx].set_priority();

					occ[q].p.unset_priority();

					if dbgf!("{:#?}: if blacklisted former occ requires no extension.", former.blacklisted()) {
						break;
					} // else need to update former.
				} else if dbgf!("{:#?}: extend(), still possible?", occ[q].extend()) {
					occ[q].mark.reset();
					// not done yet with occ[q]. It was extended and needs to be reevaluated for new mark.
					continue;
				} else {
					ks.kmp[mark_idx].blacklist();
					break;
				}

				// FIXME: can we circumvent a linear search here?
				eprintln!("search whether there was already an occurance for former position");
				q = occ.len();
				for i in 0..occ.len() {
					eprintln!("{:#x} vs {:#x}", occ[i].mark.p, former);
					if (i == 0 && occ[i].mark.p.pos() == former.pos()) || occ[i].mark.p.extpos() == former.extpos() {
						eprintln!("found {}.", i);
						q = i;
						break;
					}
					if q == occ.len() && occ[i].p.pos() == 0 && occ[i].all_kmers() {
						assert_eq!(occ[i].p.pos(), 0);
						q = i;
					}
				}
				if q == occ.len() || if q != 0 {occ[q].mark.p != former} else {occ[q].mark.p.pos() != former.pos()} {
					eprintln!("not found, we have to build occ.");
					let plimits = self.get_plimits(ks, &kc, former.b2pos());

					let is_template = (former & 1) == 0; // XXX likely incorrect. In what direction
					// do we need to extend?
					//if former.same_ori(occ[q].mark.p) {occ[q].is_template} else {!occ[q].is_template};

					let mut new_occ = Occurrence::new(plimits, kc, former.extension(), dbg!(is_template));
					eprintln!("pos: {:#x}, x:{}", former.b2pos(), former.x());
					while {
						let b2 = ks.b2_for_p(new_occ.p.pos());
						!new_occ.complete(b2) || new_occ.mark.p.b2pos() != former.b2pos()//ks, 
					} {
						debug_assert!(new_occ.p.pos() < new_occ.plim, "{:#x} < {:#x} ?", new_occ.p.pos(), new_occ.plim);
					}
					println!( "// occ insert: q:{}, {:#x}", q, new_occ.mark.p);
					if q == occ.len() {
						occ.push(new_occ);
					} else {
						occ[q] = new_occ;
					}

				} else {
					// found.
					assert!(occ[q].all_kmers(), "{}, {}", occ[q].i, kc.readlen);
					//occ[q].p.set_priority();
				}
			}
		}
		self.finalize_n_stretch_if_pending(ks);
		occ[0].p.pos() as u64
	}
}

#[cfg(test)]
mod tests {
	use super::KmerConst;
	use super::KmerIter;
	use super::KmerStore;
	const READLEN: usize = 16;
	const SEQLEN: usize = 250;
	const EXTENT: u32 = 48;

	/*#[test]
	fn test_16_c() {
		let c = KmerConst::new(READLEN, SEQLEN, extent);
		let mut ks = KmerStore::<u64>::new(c.bitlen);
		assert_eq!(ks.kmp.len(), 128);
		{
			let mut kmi = KmerIter::new(&mut ks);

			let mut vq = vec![Occurrence::new((0, u64::max_value()), &c, 0, true)];

			let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCC"[..].to_owned();
			assert_eq!(kmi.markcontig(&mut vq, &mut seq1.iter(), &c), 0x20);
			assert_eq!(vq.len(), 1);

			assert_eq!(vq[0].kmer.dna, 0x55);
			assert_eq!(vq[0].kmer.get_idx(true), 0);
		}
		assert_eq!(ks.kmp[0], (1 << c.priority_shft) | (16 << 1) - c.kmerlen as u64);
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0);
		}
		//let seq2: Vec<u8> = b"ATCGTCACTGATATCGATCC"[..].to_owned();
		//assert_eq!(kmi.markcontig(&mut vq, &mut seq2.iter(), &c), 0x48);
	}*/
	#[test]
	fn test_16n() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new();
			let seq1: Vec<u8> = b"NNNNNNNNNNNNNNNN"[..].to_owned();
			kmi.markcontig::<u64>(&mut ks, &mut seq1.iter(), &kc);
		}
		assert_eq!(ks.contig.len(), 1);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 16);
	}
	#[test]
	fn test_1n() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new();
			let seq1: Vec<u8> = b"N"[..].to_owned();
			kmi.markcontig::<u64>(&mut ks, &mut seq1.iter(), &kc);
		}
		assert_eq!(ks.contig.len(), 1);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 1);
	}
	#[test]
	fn test_1n1c1n() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		{
			let mut kmi = KmerIter::new();
			let seq1: Vec<u8> = b"NCN"[..].to_owned();
			kmi.markcontig::<u64>(&mut ks, &mut seq1.iter(), &kc);
		}
		assert_eq!(ks.contig.len(), 2);
		assert_eq!(ks.contig[0].twobit, 0);
		assert_eq!(ks.contig[0].genomic, 1);
		assert_eq!(ks.contig[1].twobit, 2);
		assert_eq!(ks.contig[1].genomic, 3);

	}
	#[test]
	fn test_17_c() {
		let kc = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::<u64>::new(kc.bitlen);
		assert_eq!(ks.kmp.len(), 128);
		{
			let mut kmi = KmerIter::new();

			// XXX ik verwacht hier (17 x C!!) een multimapper, geen positie, maar geblacklist.
			let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCC"[..].to_owned();
			kmi.markcontig::<u64>(&mut ks, &mut seq1.iter(), &kc);
		}
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		let unresolveable = (1 << 63) | ((kc.ext_max() as u64) << 48);
		assert_eq!(ks.kmp[0], unresolveable, "test_17_c(): {:x}", ks.kmp[0]);
	}
}
