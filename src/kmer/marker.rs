extern crate bio;
extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate kmerstore;
extern crate extqueue;

use self::kmerstore::KmerStore;
use self::extqueue::ExtQueue;

pub struct KmerIter<'a> {
    pos_mask: u64, // : all these are used only once.
    pos_ori_bitcount: u32, //
    priority_shft: u32, //
    kmerlen: usize,
    readlen: usize,
    pub ks: &'a mut KmerStore,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(ks: &'a mut KmerStore, readlen: usize, bitlen: usize, pos_ori_bitcount: u32) -> Self {
	KmerIter {
            pos_mask: (1 << pos_ori_bitcount) - 1,
            pos_ori_bitcount,
            priority_shft: 63,
            kmerlen: bitlen / 2,
            readlen,
            ks,
        }
    }

    // with 16 Nt we can make an entry, but minima / maxima are determined only with
    // 64 Nts. Minima for kmers with no extension; for extensions 2, 4, 8, 16, 32 and 48
    // a hashed kmer maximum is used as index, instead. (with 48 there is one entry).
    // If still not unique, include all kmers wrapping sum hashes for minima within
    // 64bp, until 1 << 16 This is to enable paired-end matching.
    pub fn markcontig(& mut self, seq: &[u8]) -> u64 {

        // position to b2 offset: Ns and contig offset excluded.
        let mut eq = ExtQueue::new(self.ks, 0, self.readlen, self.kmerlen, self.pos_ori_bitcount, true);
        //let mut eq = ExtQueue::new(0, self.readlen, self.kmerlen, true);
        let goffs = eq.loc.p & !1;
        let mut i = 0;

	for c in seq {
            // store sequence in twobit
            let mut b2 = 4;
            let p = eq.loc.p as usize >> 3;
            if let Some(qb) = eq.ks.b2.get_mut(p) {
                b2 = (*c >> 1) & 0x7; // convert ascii to 2bit
                if b2 < 4 {
                    *qb |= b2 << (eq.loc.p & 6);
                }
            }
            if b2 < 4 {
                eq.next(b2); // adds to dna/rc of kmer and increments location
                //println!("i:{}", i);
                if i >= self.kmerlen {
                    let _ = eq.progress();
                }
                i += 1;
            } else {
                eq.update_contig(goffs, self.pos_mask ^ 1);
                //eq.ks.update_contig(&mut eq, goffs, self.pos_mask ^ 1);
                i = 0;
            }
        }
        eq.loc.p
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

