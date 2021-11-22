use serde::{Deserialize, Serialize};
use std::cmp::Ordering::{Equal, Greater, Less};

use crate::kmerloc::PriExtPosOri;
//use extqueue::ExtQueue;
use anyhow::{anyhow, ensure, Result};

#[derive(Serialize, Deserialize)]
pub struct Contig {
    pub twobit: u64,
    pub genomic: u64, // if 0: next contig
}

#[derive(Serialize, Deserialize)]
pub struct KmerStore<T> {
    pub p_max: u64,
    pub opt: u64,
    pub b2: Vec<u8>,
    pub kmp: Vec<T>, // position + strand per k-mer.
    pub contig: Vec<Contig>,
}

impl<T: PriExtPosOri> KmerStore<T> {
    pub fn new(bitlen: usize) -> Self {
        let shift = bitlen - 2;
        KmerStore {
            p_max: 0,
            opt: 0,
            b2: vec![0; 1 << shift], // sequence (4 per u8).
            kmp: vec![T::from_u64(0).unwrap(); 1 << (shift + 1)], // kmer positions
            contig: Vec::new(),      // contig info
        }
    }
    pub fn push_contig(&mut self, p: u64, goffs: u64) {
        self.contig.push(Contig {
            twobit: p, // TODO: bits could mean something about contig.
            genomic: goffs,
        });
    }
    /// Adjust offset for contig. The
    pub fn offset_contig(&mut self, offset: u64) {
        if let Some(ctg) = self.contig.last_mut() {
            ctg.genomic += offset;
        }
    }
    pub fn get_contig(&self, p: u64) -> usize {
        let mut size = self.contig.len();
        let mut base = 0;
        while size > 0 {
            size /= 2;
            let mid = base + size;
            base = match (self.contig[mid].twobit).cmp(&p) {
                Less => mid,
                Greater | Equal => base,
            };
        }
        base
    }
    pub fn get_twobit_after(&self, i: usize) -> Result<u64> {
        ensure!(
            i + 1 < self.contig.len(),
            "get_twobit_after(): End of contig"
        );
        Ok(self.contig[i + 1].twobit)
    }
    pub fn get_twobit_before(&self, i: usize) -> Result<u64> {
        ensure!(i != 0, "get_twobit_before(): Start of contig");
        Ok(self.contig[i - 1].twobit)
    }
    pub fn b2_for_p(&self, p: u64) -> Result<u8> {
        let i = p as usize >> 3;
        self.b2
            .get(i)
            .map(|x| (x >> (p & 6)) & 3)
            .ok_or_else(|| anyhow!("b2_for_p():OOB {:x} {:x}", p, self.b2.len()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rdbg::STAT_DB;
    #[test]
    fn it_works() {
        //let mut ks = KmerStore::new(32); // allocates 16 gig of data
        let mut ks = KmerStore::<u64>::new(2);
        let mut p = 0;
        let goffs = p & !1;

        ks.push_contig(p, goffs);
        ks.offset_contig(10000); // simulate N-stretch of 10000
        p += 64; // 64 Nt's == one readlength
        ks.push_contig(p, goffs);
        ks.offset_contig(64);

        ks.push_contig(p, goffs);
        ks.offset_contig(10000); // another N-stretch of 10000
        ks.push_contig(p, goffs);
        let i = ks.get_contig(11016);
        dbg_assert_eq!(ks.contig[i].twobit, 64);
    }
}
