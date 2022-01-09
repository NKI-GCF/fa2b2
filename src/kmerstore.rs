use crate::kmer::{TwoBit, TwoBitx4};
use crate::kmerloc::ExtPosEtc;
use crate::new_types::position::Position;
use crate::rdbg::STAT_DB;
use ahash::AHashMap;
use anyhow::{anyhow, ensure, Result};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering::{Equal, Greater, Less};

#[derive(Serialize, Deserialize)]
pub struct Contig {
    pub twobit: Position, // position required for binary search
    pub genomic: u64,     // if 0: next contig
}

type ContigRng = (Position, Position);
pub type Repeat = (u32, u32);

#[derive(Serialize, Deserialize)]
pub struct KmerStore {
    pub pos_max: Position,
    pub opt: u64,
    pub rep_max_dist: Position,
    pub b2: Vec<u8>,
    pub kmp: Vec<u64>, // position + strand per k-mer.
    pub contig: Vec<Contig>,
    pub repeat: AHashMap<Position, Repeat>,
}

impl KmerStore {
    pub fn new(bitlen: usize, rep_max_dist: u64) -> Self {
        let shift = bitlen - 2;
        KmerStore {
            pos_max: Position::zero(),
            opt: 0,
            rep_max_dist: rep_max_dist.basepos_to_pos(),
            b2: vec![0; 1 << shift],                  // sequence (4 per u8).
            kmp: vec![u64::zero(); 1 << (shift + 1)], // kmer positions
            //kmp,
            contig: Vec::new(), // contig info
            repeat: AHashMap::new(),
        }
    }
    pub fn store_twobit(&mut self, pos: Position, b2: TwoBit) {
        self.b2[pos.byte_pos()] |= b2.pos_shift(pos).as_u8();
    }
    pub fn push_contig(&mut self, pos: Position, goffs: u64) {
        self.contig.push(Contig {
            twobit: pos,
            genomic: goffs,
        });
    }
    /// Adjust offset for contig.
    pub fn offset_contig(&mut self, offset: u64) {
        if let Some(ctg) = self.contig.last_mut() {
            ctg.genomic += offset;
        }
    }
    pub fn set_kmp(&mut self, min_idx: usize, min_p: ExtPosEtc) {
        self.kmp[min_idx].set(min_p);
    }
    pub fn get_bitlen(&self) -> usize {
        self.b2.len().trailing_zeros() as usize + 2
    }

    /// binary search contig lower boundary
    fn get_contig(&self, pos: Position) -> usize {
        let mut size = self.contig.len();
        let mut base = 0;
        while size > 0 {
            size /= 2;
            let mid = base + size;
            base = match (self.contig[mid].twobit).cmp(&pos) {
                Less => mid,
                Greater | Equal => base,
            };
        }
        base
    }
    pub fn get_contig_start_end_for_p(&self, pos: Position) -> ContigRng {
        let i = self.get_contig(pos);

        dbg_assert!(i < self.contig.len());
        dbg_assert!(self.contig[i].twobit <= pos);
        (
            self.get_twobit_before(i).unwrap_or(Position::zero()),
            self.get_twobit_after(i).unwrap_or(self.pos_max),
        )
    }
    fn get_twobit_after(&self, i: usize) -> Result<Position> {
        ensure!(
            i + 1 < self.contig.len(),
            "get_twobit_after(): End of contig"
        );
        Ok(self.contig[i + 1].twobit)
    }
    fn get_twobit_before(&self, i: usize) -> Result<Position> {
        ensure!(i != 0, "get_twobit_before(): Start of contig");
        Ok(self.contig[i - 1].twobit)
    }
    pub fn b2_for_p(&self, p: ExtPosEtc, is_repeat: bool) -> Result<TwoBit> {
        let pos = p.pos();
        ensure!(
            pos < self.pos_max || is_repeat,
            "running into sequence head"
        );
        self.b2
            .get(pos.byte_pos())
            .map(|b2x4| TwoBitx4::from(b2x4).to_b2(p, is_repeat))
            .ok_or_else(|| anyhow!("stored pos past contig? {:#x}", p))
    }
    pub fn extend_repetitive(&mut self, min_pos: Position, dist: u64) {
        let dist_u32 = u32::try_from(dist).expect("dist for repeat extension doesn't fit in u32");
        let repeat = self.repeat.entry(min_pos).or_insert((dist_u32, 0));
        repeat.1 = dist_u32;
    }

    pub fn replace_repetitive(&mut self, old_pos: Position, new_pos: Position, dist: u64) {
        let dist_u32 = u32::try_from(dist).expect("dist for repeat extension doesn't fit in u32");
        let old_repeat = self.repeat.remove(&old_pos).expect("oldpos not stored");

        let repeat = self.repeat.entry(new_pos).or_insert(old_repeat);
        repeat.1 += dist_u32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rdbg::STAT_DB;
    #[test]
    fn it_works() {
        //let mut ks = KmerStore::new(32); // allocates 16 gig of data
        let mut ks = KmerStore::new(2, 10_000);
        let mut p = 0;
        let goffs = 0;

        ks.push_contig(p.pos(), goffs);
        ks.offset_contig(10_000); // simulate N-stretch of 10000
        p += 64.basepos_to_pos().as_u64(); // one readlength
        ks.push_contig(p.pos(), goffs);
        ks.offset_contig(64);

        ks.push_contig(p.pos(), goffs);
        ks.offset_contig(10_000); // another N-stretch of 10000
        ks.push_contig(p.pos(), goffs);

        let i = ks.get_contig(5_508.basepos_to_pos());
        dbg_assert_eq!(ks.contig[i].twobit, 64.basepos_to_pos());
    }
}
