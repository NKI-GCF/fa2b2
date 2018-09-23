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
            let mut eq = &mut vq[y];
            let b2;
            if y == 0 {
                match seq.next() {
                    None => break, //                                                    loop exit loop
                    Some(c) => {
                        b2 = (*c >> 1) & 0x7;
                        if b2 >= 4 {
                            eq.update_non_contig(self.ks, goffs);
                            continue;
                        }
                        self.add_new_seq(eq.loc.p, b2); //                                 store twobit
                    },
                }
            } else {
                b2 = self.ks.b2_for_p(eq.loc.p);
            }
            if !eq.next(b2) { //                                         false if kmer not yet complete
                continue;
            }
            while eq.x <= eq.ext {
                let new = eq.eval_hash_turnover();
                let ext_bits = (eq.x << kc.pos_ori_bitcount) as u64 | eq.priority;

                let former = self.ks.kmp[new.idx];
                if !eq.is_kmer_available(new.p, former, ext_bits) {

                    let in_use = former & (1 << kc.priority_shft);
                    let n = (former ^ in_use) >> kc.pos_ori_bitcount;

                    if eq.extension(former) > ext_bits { //                  stored has higher priority
                        eq.x += 1;
                        continue;
                    }
                    self.ks.kmp[new.idx] = in_use | ((n + 1) << kc.pos_ori_bitcount);

                    // XXX: plimits end is not restrictive enough at end of seq!
                    plimits = eq.get_plimits(self.ks, eq.b2pos(former));
                    is_forward = eq.is_rebuild_forward(former, new.p);

                    println!("[{}, {}, {}]: {:x} in use by {:x} {:?}, {}", eq.i, eq.x, y,
                             new.idx, former, plimits, is_forward);

                    // FIXME: the position may already be in vq. then rebuild is not required.
                    // but we need to temporarily switch to that extension

                    y += 1;
                    println!("look back no. {}", y);
                    assert!(y < 1000);
                    break;
                }
                self.ks.kmp[new.idx] = new.p | ext_bits;
                println!("[{}, {}, {}]: {:x} set to 0x{:x}", eq.i, eq.x, y, new.idx, new.p | ext_bits);
                eq.priority = 0;
                if y != 0 && eq.all_kmers() {
                    y -= 1;
                    break;
                }
                // if there is one recurring, there may be more, therefore we should
                // store the entire eq. TODO: eq by eq extension: problem: orientation.
                // rev_relocate, afh van orientatie van oldp & 1?
                eq.marked[eq.x as usize].push_back(new);
                eq.x += 1;
            }
            if eq.priority != 0 {
                // TODO:
                // If still not unique, include all kmers wrapping sum hashes for minima within
                // 64bp, until 1 << 16; to enable paired-end matching.
                //
                // Maybe, not every, but only for neighbouring unmappable sites?
            }
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

