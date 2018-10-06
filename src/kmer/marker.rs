extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate occurrence;
extern crate kmerconst;
extern crate kmerloc;
extern crate kmerstore;

use self::occurrence::Occurrence;
use self::kmerconst::KmerConst;
use self::kmerstore::KmerStore;
use std::cmp::Ordering::{Equal, Greater, Less};
use std::slice::Iter;

pub struct KmerIter<'a> {
	pub ks: &'a mut KmerStore,
} //^-^\\

impl<'a> KmerIter<'a> {
	pub fn new(ks: &'a mut KmerStore) -> Self {
		KmerIter { ks }
	}

	/// add ascii sequence converted to 2bit to kmers and store.
	fn add_new_seq(&mut self, p: u64, b2: u8) {
		if let Some(qb) = self.ks.b2.get_mut(p as usize >> 3) {
			if b2 < 4 {
				*qb |= b2 << (p & 6);
			}
		}
		self.ks.p_max = (p + 2) & !1;
	}
	fn get_vq(&self, vq: &Vec<Occurrence<'a>>, p: u64) -> usize {
		let mut size = vq.len();
		let mut base = 0;
		while size > 0 {
			size /= 2;
			let mid = base + size;
			base = match (vq[mid].loc.p).cmp(&p) {
				Less => mid,
				Greater | Equal => base,
			};
		}
		base
	}

	fn next(&mut self, seq: &mut Iter<u8>, occ: &mut Occurrence<'a>, goffs: u64) -> Option<u8> {
		seq.next().and_then(|c| Some((*c >> 1) & 0x7)).map(|b2| {
			if b2 >= 4 {
				occ.update_non_contig(self.ks, goffs);
			} else {
				self.add_new_seq(occ.loc.p, b2); // store twobit
				occ.x = 0;
			}
			b2
		})
	}

	/// Traverse a new reference sequence. Store as twobit, and update extension queue,
	/// N-stretches excluded; Contig regions indicate chrosomes start and the offset of such N-stretches.
	pub fn markcontig(
		&mut self,
		occ: &mut Vec<Occurrence<'a>>,
		seq: &mut Iter<u8>,
		kc: &'a KmerConst,
	) -> u64 {
		// position to b2 offset: Ns and contig offset excluded.
		let mut y = 0;
		let mut lastp: Vec<u64> = vec![];
		let mut plimits = (occ[y].loc.p, occ[y].plim);
		let mut is_forward = occ[y].is_template;
		let mut extend_kmer = false;
		let mut former = 0;
		let goffs = plimits.0 & !1;
		occ[y].update_non_contig(self.ks, goffs);
		// vq kan gesorteerd op loc.p, en dan binary gesearched worden.

		loop {
			if !extend_kmer {
				if y != 0 && occ[y].loc.p >= self.ks.p_max {
					panic!("need more seq (returning to lookback no. 0)!");
					//y = 0;
				}
				let b2 = if y == 0 {
					match self.next(seq, &mut occ[0], goffs) {
						Some(7) => continue,
						Some(b2) => b2,
						None => break,
					}
				} else {
					self.ks.b2_for_p(occ[y].loc.p)
				};
				if !occ[y].next(b2) {
					//false if kmer not yet complete
					continue;
				}
			} else {
				extend_kmer = false;
			}
			let mut x = occ[y].x;
			while (1 << x) <= occ[y].ext {
				let new = occ[y].eval_hash_turnover();

				println!("ixyp:[{}, {}, {}, {:x}]", occ[y].i, x, y, occ[y].loc.p);
				let ext_bits = (x << kc.extent) as u64 | occ[y].priority;

				former = self.ks.kmp[new.idx];
				if !occ[y].is_kmer_available(new.p, former, ext_bits) {
					let in_use = former & (1 << kc.priority_shft);
					let n = (former ^ in_use) & !kc.pos_mask;
					x += 1;
					occ[y].x = x;

					if occ[y].extension(former) > ext_bits {
						//                  stored has higher priority
						println!(
							"stored has higher priority (or blacklisted), {:x}, {:x}",
							former, ext_bits
						);
						continue;
					}
					lastp.push(occ[y].loc.p);
					println!(
						"blacklisting kmer, vq[{}].priority == {:x}",
						y, occ[y].priority
					);
					extend_kmer = true;
					self.ks.kmp[new.idx] = in_use | (n + (1 << kc.extent)); //         blacklist

					plimits = occ[y].get_plimits(self.ks, occ[y].b2pos(former));
					is_forward = occ[y].is_rebuild_forward(former, new.p);
					break;
				}
				if (former & kc.pos_mask) == 0 {
					self.ks.kmp[new.idx] = (new.p - kc.kmerlen as u64) | ext_bits;
					occ[y].priority = 0;
					former = lastp.pop().unwrap_or(0);
					//y = if let Some(p) = lastp.pop() {self.get_vq(&vq, p)} else {0};
					//println!("set; y: {}, vq[y].priority == {:x}", y, vq[y].priority);
					break;
				} else if !occ[y].is_within_read(new.p | ext_bits, former) {
					//panic!("check branch");
					if occ[y].extension(former) == ext_bits {
						extend_kmer = true;
						println!(
							"equal extension: switch to last extension, without pop {}, {}",
							x,
							occ[y - 1].x
						);
						while {
							assert_ne!(y, 0);
							y -= 1;
							occ[y].loc.p != occ[y].b2pos(former)
						} {}
						x += 1;
						occ[y].x = x;
						break;
					}
					self.ks.kmp[new.idx] = (new.p - kc.kmerlen as u64) | ext_bits;
					occ[y].priority = 0;
					if y != 0 {
						extend_kmer = true;
					}
					println!(
						"[{}, {}, {}, {:x}]: {:x} set to 0x{:x}",
						occ[y].i, x, y, occ[y].loc.p, new.idx, self.ks.kmp[new.idx]
					);
				}
				if occ[y].priority == 0 && (self.ks.opt & 2) == 0 {
					println!("leaving {}", y);
					break;
				}
				// if there is one recurring, there may be more, therefore we should
				// store the entire eq. TODO: eq by eq extension: problem: orientation.
				// rev_relocate, afh van orientatie van oldp & 1?
				occ[y].marked[x as usize].push_back(new);
				x += 1;
				occ[y].x = x;
			}
			if (former & kc.pos_mask) != 0 {
				extend_kmer = true;
				y = self.get_vq(&occ, former);
				if occ[y].loc.p != former {
					println!("vq.insert({}, Occurrence[{:?}, {}])", y, plimits, x);
					occ.insert(y, Occurrence::new(plimits, kc, x, is_forward));
					y += 1;
					former = 0;
				}
			}
			if occ[y].ext > occ[y].ext_max && occ[y].priority != 0 {
				println!("not able to complete lookback: {}, p: {:x}", y, occ[y].loc.p);
				// TODO:
				// If still not unique, include all kmers wrapping sum hashes for minima within
				// 64bp, until 1 << 16; to enable paired-end matching.
				//
				// Maybe, not every, but only for neighbouring unmappable sites?
			}
		}
		occ[y].loc.p & !1
	}
}

#[cfg(test)]
mod tests {
	use super::Occurrence;
	use super::KmerConst;
	use super::KmerIter;
	use super::KmerStore;
	const READLEN: usize = 16;
	const SEQLEN: usize = 250;
	const extent: u32 = 48;

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
		let c = KmerConst::new(READLEN, SEQLEN, extent);
		let mut ks = KmerStore::new(c.bitlen);
		assert_eq!(ks.kmp.len(), 128);
		{
			let mut kmi = KmerIter::new(&mut ks);

			let mut vq = vec![Occurrence::new((0, u64::max_value()), &c, 0, true)];

			// XXX ik verwacht hier (17 x C!!) een multimapper, geen positie, maar geblacklist.
			let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCC"[..].to_owned();
			//let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..].to_owned();
			//let seq1: Vec<u8> = b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..].to_owned();
			kmi.markcontig(&mut vq, &mut seq1.iter(), &c);
			if let Some(eq) = vq.get(0) {
				assert_eq!(eq.kmer.dna, 0x55); // index out of bounds????
				assert_eq!(eq.kmer.get_idx(true), 0);
			} else {
				panic!("out of bounds!?");
			}
		}
		//800100000000001a
		assert_ne!(ks.kmp[0], 0);
		for i in 1..ks.kmp.len() {
			assert_eq!(ks.kmp[i], 0, "[{}], {:x}", i, ks.kmp[i]);
		}
		//let seq2: Vec<u8> = b"ATCGTCACTGATATCGATCC"[..].to_owned();
		//assert_eq!(kmi.markcontig(&mut vq, &mut seq2.iter(), &c), 0x48);
	}
}
