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

pub struct KmerIter<'a> {
	pub ks: &'a mut KmerStore,
	n_stretch: u64,
} //^-^\\

impl<'a> KmerIter<'a> {
	pub fn new(ks: &'a mut KmerStore) -> Self {
		KmerIter { ks, n_stretch: 0 }
	}

	fn complete_occurance(&mut self, occ: &mut Occurrence, b2: u8) -> bool {
		let p = occ.loc.p.low();
		if let Some(qb) = self.ks.b2.get_mut(p as usize >> 3) {
			*qb |= b2 << (p & 6);
		}
		self.ks.p_max = (p + 2) & !1;
		occ.add(b2);
		if occ.all_kmers() {
			occ.turnover();
			true
		} else if occ.now_all_kmers() {
			self.ks
				.offset_contig(self.n_stretch + 1 - occ.kc.readlen as u64);
			true
		} else {
			false
		}
	}

	/// if N, insert contig (once) and count stretch, false
	/// else store and update occurrance kmers, complete ? true : false
	fn complete_occurance_or_contig(
		&mut self,
		occ: &mut Occurrence<'a>,
		b2: u8,
		goffs: u64,
	) -> bool {
		if b2 < 4 && self.complete_occurance(occ, b2) {
			true
		} else {
			// N-stretch is incremented even until kmer is finished, this is corrected when incrementing contig.
			if self.n_stretch == 0 {
				let p = occ.loc.p.low() - 2; //was already incremented in complete_occurance() -> occ.add() -> loc.next()
				self.ks.push_contig(p, goffs);

				// clear all except orientation and position to rebuild at the start of a new contig.
				assert_eq!(occ.i, 1);
				occ.clear();
			}
			self.n_stretch += 1;
			false
		}
	}

	/*fn build_occ(&mut self, kc: &KmerConst, former: &KmerLoc, x: usize) -> Occurrence {
		let plimits = self.get_plimits(&kc, former.b2pos());

		println!( "// not found. {:016x} plim:{:?} {}", former.get(), plimits, x);

		let is_forward = true; // FIXME
		let mut new_occ = Occurrence::new(plimits, kc, x as u64, is_forward);
		while {
			let b2 = self.ks.b2_for_p(new_occ.loc.low());
			!self.complete_occurance(&mut new_occ, b2)
				|| new_occ.marked[x][0].b2pos() != former.b2pos()
		} {}
		println!("// occ insert: {:x}", new_occ.get_extreme().get());
		new_occ
	}*/

	fn get_occ_idx<T: MidPos>(&self, occ: &Vec<Occurrence<'a>>, former: &T) -> usize {

		println!("search");
		let mut first_unused = occ.len();
		for i in 0..occ.len() {
			println!("({:x}){:x} vs {:x}", occ[i].loc.p, occ[i].get_extreme().p, former.to_u64().unwrap());
			if occ[i].get_extreme().p.extpos() == former.extpos() {
				println!("found.");
				return i;
			}
			if first_unused == occ.len()
				&& occ[i].loc.p.top() == occ[i].ext
				&& occ[i].loc.priority() != 0
			{
				assert_eq!(occ[i].loc.p.low(), 0);
				first_unused = i;
			}
		}
		if first_unused == occ.len() {
			println!("not found.");
		}
		return first_unused;
	}

	/// TODO: replace this function
	/// when rebuilding eq for recurrent kmer, and extending take into account contig boundaries
	/// for that site
	pub fn get_plimits(&self, kc: &'a KmerConst, p: u64) -> (u64, u64) {
		// binary search; limit endp to end of contig
		assert!(p != 0);
		let i = self.ks.get_contig(p);
		assert!(i < self.ks.contig.len());
		//println!("contig[{}].twobit:{}, p:{}", i, ks.contig[i].twobit, p);

		assert!(self.ks.contig[i].twobit <= p);
		let x = self.ks.get_twobit_before(i).unwrap_or(0);
		let addl = (kc.max_no_kmers << 1) as u64;
		let p_max = self.ks.p_max;
		let left = if p > x + addl { p - addl } else { x };
		let right = cmp::min(
			p + (kc.readlen << 1) as u64,
			self.ks.get_twobit_after(i).unwrap_or(p_max),
		);
		(left, right)
	}
	pub fn markcontig<T: MidPos>(
		&mut self,
		occ: &mut Vec<Occurrence<'a>>,
		seq: &mut Iter<u8>,
		kc: &'a KmerConst,
	) -> u64 {
		let mut q = 0;
		let goffs = occ[0].loc.p.low();
		let mut former_stack: Vec<u64> = vec![];
		let mut previous = None;

		while let Some(c) = seq.next() {
			if !self.complete_occurance_or_contig(&mut occ[0], (*c >> 1) & 0x7, goffs) {
				continue;
			}
			while occ[q].loc.p.extension() <= occ[q].ext {

				let new = occ[q].eval_hash_turnover();
				if previous == Some(new) {
					println!("same: [{}][{:08x}] = {:016x}", q, new.idx, self.ks.kmp[new.idx]);
					//break; or continue?
				}
				previous = Some(new);

				let ext_bits = occ[q].loc.p.top();

				let mut former: T = T::from_u64(self.ks.kmp[new.idx]).unwrap();

				if former.to_u64().unwrap() <= ext_bits {
					self.ks.kmp[new.idx] = new.p | (1 << 63);
					println!( "// [{}][{:08x}] = {:016x}", q, new.idx, self.ks.kmp[new.idx]);

					occ[q].loc.unset_priority();
					if former.b2pos() == 0 {
						former = match former_stack.pop() {
							Some(f) => T::from_u64(f).unwrap(),
							None => break,
						};
						println!("\n// retrieved {:x} from former_stack", former.to_u64().unwrap());
					}
				} else {
					occ[q].loc.p.extend();
					if occ[q].loc.p.extension() > occ[q].ext {
						println!("unresolveable: TODO: kmers of adjoining regions like these?");
						break;
					}

					if former.b2pos() == 0 || former.top() > ext_bits {
						println!( "// kept former: blacklist for or entry with higher extension.");
						previous = None;
						continue;
					}
					// identical extension: blacklist, handle former and then this p
					self.ks.kmp[new.idx] += 1 << kc.extent;
					self.ks.kmp[new.idx] &= !kc.pos_mask;
					println!("// unreplaceable, former_stack({}) push: occ[{}].loc.p:{:016x}", former_stack.len(), q, occ[q].loc.p);
					former_stack.push(occ[q].get_extreme().p);
				}
				let oldp = occ[q].loc.p.low();
				q = self.get_occ_idx::<T>(&occ, &former);
				if q == occ.len() || occ[q].get_extreme().p != former.to_u64().unwrap() {
					let plimits = self.get_plimits(&kc, former.b2pos());
					let x = former.x();

					// not found.

					let is_forward = true; // FIXME
					let mut new_occ = Occurrence::new(plimits, kc, x as u64, is_forward);
					while {
						let b2 = self.ks.b2_for_p(new_occ.loc.p.low());
						!self.complete_occurance(&mut new_occ, b2)
							|| new_occ.marked[x][0].p.b2pos() != former.b2pos()
					} {}
					println!( "// occ insert: q:{}, {:x}", q, new_occ.get_extreme().p);
					if q == occ.len() {
						occ.push(new_occ);
					} else {
						occ[q] = new_occ;
					}

				} else {
					// found.
					occ[q].loc.set_priority();
				}
			}
		}
		occ[0].loc.p.low() - kc.kmerlen as u64
	}
}

#[cfg(test)]
mod tests {
	use super::KmerConst;
	use super::KmerIter;
	use super::KmerStore;
	use super::Occurrence;
	const READLEN: usize = 16;
	const SEQLEN: usize = 250;
	const EXTENT: u32 = 48;

	/*#[test]
    fn test_16_c() {
        let c = KmerConst::new(READLEN, SEQLEN, extent);
        let mut ks = KmerStore::new(c.bitlen);
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
	fn test_17_c() {
		let c = KmerConst::new(READLEN, SEQLEN, EXTENT);
		let mut ks = KmerStore::new(c.bitlen);
		assert_eq!(ks.kmp.len(), 128);
		{
			let mut kmi = KmerIter::new(&mut ks);

			let mut vq = vec![Occurrence::new((0, u64::max_value()), &c, 0, true)];

			// XXX ik verwacht hier (17 x C!!) een multimapper, geen positie, maar geblacklist.
			let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCC"[..].to_owned();
			//let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..].to_owned();
			//let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..].to_owned();
			kmi.markcontig::<u64>(&mut vq, &mut seq1.iter(), &c);
			if let Some(eq) = vq.get(0) {
				assert_eq!(eq.kmer.dna, 0x55); // index out of bounds????
				assert_eq!(eq.kmer.get_idx(true), 0);
			} else {
				panic!("out of bounds!?");
			}
		}
		//800100000000001a
		println!("{:x}", ks.kmp[0]);
		assert_ne!(ks.kmp[0], 0);
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		//let seq2: Vec<u8> = b"ATCGTCACTGATATCGATCC"[..].to_owned();
		//assert_eq!(kmi.markcontig(&mut vq, &mut seq2.iter(), &c), 0x48);
	}
}
