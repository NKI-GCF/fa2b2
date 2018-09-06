extern crate bio;
extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate kmerconst;
extern crate extqueue;

use self::extqueue::ExtQueue;
use self::kmerconst::KmerConst;


// with 16 Nt we can make an entry, but minima / maxima are determined only with
// 64 Nts. Minima for kmers with no extension; for extensions 2, 4, 8, 16, 32 and 48
// a hashed kmer maximum is used as index, instead. (with 48 there is one entry).
// If still not unique, include all kmers wrapping sum hashes for minima within
// 64bp, until 1 << 16 This is to enable paired-end matching.
pub fn markcontig(vq: &mut Vec<ExtQueue>, seq: &[u8], kc: &KmerConst) -> u64 {

    // position to b2 offset: Ns and contig offset excluded.
    //let mut eq = ExtQueue::new(self.ks, 0, self.readlen, self.kmerlen, self.pos_ori_bitcount, true);
    let goffs = vq[0].loc.p & !1;
    let mut i = 0;

    for c in seq {
        // store sequence in twobit
        let mut b2 = 4;
        let p = vq[0].loc.p as usize;
        if let Some(qb) = vq[0].ks.b2.get_mut(p >> 3) {
            b2 = (*c >> 1) & 0x7; // convert ascii to 2bit
            if b2 < 4 {
                *qb |= b2 << (p & 6);
            }
        }
        if b2 < 4 {
            vq[0].next(b2); // adds to dna/rc of kmer and increments location
            //println!("i:{}", i);
            if i >= kc.kmerlen {
                let _ = vq[0].progress();
            }
            i += 1;
        } else {
            vq[0].update_contig(goffs);
            //vq[0].ks.update_contig(&mut vq[0], goffs, pos_mask ^ 1);
            i = 0;
        }
    }
    vq[0].loc.p
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

