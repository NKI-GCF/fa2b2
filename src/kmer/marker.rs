extern crate bio;
extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate kmerconst;
extern crate kmerstore;
extern crate extqueue;
extern crate kmerloc;

use self::extqueue::ExtQueue;
use self::kmerconst::KmerConst;
use self::kmerloc::KmerLoc;
use self::kmerstore::KmerStore;
use std::slice::Iter;

pub struct KmerIter<'a> {
    pub ks: &'a mut KmerStore,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore) -> Self {
	KmerIter {
            ks,
        }
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

    /// Traverse a new reference sequence. Store as twobit, and update extension queue,
    /// N-stretches excluded; Contig regions indicate chrosomes start and the offset of such N-stretches.
    pub fn markcontig(&mut self, vq: &mut Vec<ExtQueue<'a>>, seq: &mut Iter<u8>, kc: &'a KmerConst) -> u64 {

        // position to b2 offset: Ns and contig offset excluded.
        let mut y = 0;
        let mut plimits = (vq[y].loc.p, vq[y].plim);
        let mut is_forward = vq[y].is_template;
        let goffs = plimits.0 & !1;

        loop {
            if y + 2 == vq.len() {
                vq.pop();
            } else if y == vq.len() {
                vq.push(ExtQueue::new(plimits, kc, is_forward));
            }
            let b2;
            if y != 0 && vq[y].loc.p >= self.ks.p_max {
                //println!("need more seq (returning to lookback no. 0)!");
                y = 0;
            }
            if y == 0 {
                match seq.next() {
                    None => break, //                                                    loop exit loop
                    Some(c) => {
                        b2 = (*c >> 1) & 0x7;
                        if b2 >= 4 {
                            vq[y].update_non_contig(self.ks, goffs);
                            continue;
                        }
                        self.add_new_seq(vq[y].loc.p, b2); //                                 store twobit
                    },
                }
            } else {
                b2 = self.ks.b2_for_p(vq[y].loc.p);
            }
            if !vq[y].next(b2) { //                                         false if kmer not yet complete
                continue;
            }
            while vq[y].x <= vq[y].ext {
                let new = vq[y].eval_hash_turnover();
                let ext_bits = (vq[y].x << kc.pos_ori_bitcount) as u64 | vq[y].priority;

                let former = self.ks.kmp[new.idx];
                if !vq[y].is_kmer_available(new.p, former, ext_bits) {

                    let in_use = former & (1 << kc.priority_shft);
                    let n = (former ^ in_use) >> kc.pos_ori_bitcount;

                    if vq[y].extension(former) > ext_bits { //                  stored has higher priority
                        vq[y].x += 1;
                        continue;
                    }

                    self.ks.kmp[new.idx] = in_use | ((n + 1) << kc.pos_ori_bitcount);

                    // XXX: plimits end is not restrictive enough at end of seq!
                    plimits = vq[y].get_plimits(self.ks, vq[y].b2pos(former));
                    is_forward = vq[y].is_rebuild_forward(former, new.p);

                    if y == 0 && vq.len() != 1 {
                        // returned to 0 for more seq, now we may want to return to
                        // last y, if the recurrance is a continuation.
                        y = vq.len() - 1;
                        //println!("[recurrent] returning to lookback no. {}", y);
                        if vq[y].is_too_nearby(vq[y].loc.p, former) {
                            break;
                        }
                    }

                    //println!("[{}, {}, {}]: {:x} in use by {:x} {:?}, {}", vq[y].i, vq[y].x, y,
                    //         new.idx, former, plimits, is_forward);

                    // FIXME: the position may already be in vq. then rebuild is not required.
                    // but we need to temporarily switch to that extension

                    y += 1;
                    //println!("look back no. {}", y);
                    assert!(y < 1000);
                    break;
                }
                if former == 0 || !vq[y].is_too_nearby(new.p | ext_bits, former) {
                    self.ks.kmp[new.idx] = new.p | ext_bits;
                    //println!("[{}, {}, {}]: {:x} set to 0x{:x}", vq[y].i, vq[y].x, y, new.idx, new.p | ext_bits);
                }
                vq[y].priority = 0;
                if y != 0 && vq[y].all_kmers() {
                    y -= 1;
                    break;
                }
                // if there is one recurring, there may be more, therefore we should
                // store the entire eq. TODO: eq by eq extension: problem: orientation.
                // rev_relocate, afh van orientatie van oldp & 1?
                let x = vq[y].x as usize;
                vq[y].marked[x].push_back(new);
                vq[y].x += 1;
                if y == 0 && vq.len() != 1 {
                    y = vq.len() - 1;
                    //println!("[new] returning to lookback no. {}", y);
                    break;
                }
            }
            //if vq[y].priority != 0 {
                // TODO:
                // If still not unique, include all kmers wrapping sum hashes for minima within
                // 64bp, until 1 << 16; to enable paired-end matching.
                //
                // Maybe, not every, but only for neighbouring unmappable sites?
            //}
        }
        vq[y].loc.p & !1
    }
}


#[cfg(test)]
mod tests {
    use super::ExtQueue;
    use super::KmerIter;
    const READLEN: usize = 16;
    const KMERLEN: usize = 4;
    const BITLEN: usize = KMERLEN * 2;

    #[test]
    fn test_kmi() {
        let mut kmi = KmerIter::new(READLEN, BITLEN, 48);
        let seq: Vec<u8> = b"CCCCCC"[..].to_owned();
        let p = kmi.markcontig(&b"CCCCCC"[..].to_owned());
        assert_eq!(p, 0xc);
    }
}

