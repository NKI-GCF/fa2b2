extern crate bio;
extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate kmerconst;
extern crate extqueue;

use self::extqueue::ExtQueue;
use self::kmerconst::KmerConst;

pub struct KmerIter<'a> {
    vq: &'a mut Vec<ExtQueue<'a>>,
} //^-^\\

impl<'a> KmerIter<'a> {
    pub fn new(vq: &'a mut Vec<ExtQueue<'a>>) -> Self {
	KmerIter {
            vq,
        }
    }

    /// Traverse a new reference sequence. Store as twobit, and update extension queue,
    /// N-stretches excluded; Contig regions indicate chrosomes start and the offset of such N-stretches.
    pub fn markcontig(& mut self, seq: &[u8], kc: &KmerConst) -> u64 {
        let mut p = 0;

        // position to b2 offset: Ns and contig offset excluded.
        if let Some(eq) = self.vq.get_mut(0) {
            let goffs = eq.loc.p & !1;

            for c in seq {
                // store sequence in twobit
                let mut b2 = 4;
                p = eq.loc.p;
                if let Some(qb) = eq.ks.b2.get_mut(p as usize >> 3) {
                    b2 = (*c >> 1) & 0x7; // convert ascii to 2bit
                    if b2 < 4 {
                        *qb |= b2 << (p & 6);
                    }
                }
                if b2 < 4 {
                    eq.next(b2); // adds to dna/rc of kmer and increments location
                    if eq.d.len() < kc.kmerlen {
                        let _ = eq.progress();
                    }
                } else {
                    eq.update_non_contig(goffs);
                }
            }
        }
        p
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

