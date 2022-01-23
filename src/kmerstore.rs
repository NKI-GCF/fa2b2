use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{
    position::{BasePos, PosRange, Position},
    twobit::{TwoBit, TwoBitx4},
};
use crate::rdbg::STAT_DB;
use crate::xmer_location::XmerLoc;
use ahash::AHashMap;
use anyhow::{bail, Result};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct Contig {
    pub twobit: Position, // position required for binary search
    pub genomic: BasePos, // if 0: next contig
}

pub type Repeat = (u32, u32);

// XXX if repeats are regular and gt readlen, why store them in ks.b2?

#[derive(Serialize, Deserialize)]
pub struct KmerStore {
    b6pos: usize,
    pub(crate) opt: u64,
    pub(crate) b2: Vec<u8>,
    pub(crate) kmp: Vec<ExtPosEtc>, // position + strand per k-mer.
    pub(crate) contig: Vec<Contig>,
    pub(crate) repeat: AHashMap<Position, Repeat>,
    pub(crate) rep_max_dist: Position,
    pub(crate) seed: u16,
}

impl KmerStore {
    pub(crate) fn new(bitlen: usize, rep_max_dist: u32, seed: u16) -> Result<Self> {
        let shift = bitlen - 2;
        Ok(KmerStore {
            b6pos: 0,
            opt: 0,
            rep_max_dist: Position::from_basepos(rep_max_dist),
            b2: vec![0; 1 << shift], // sequence (4 per u8).
            kmp: vec![ExtPosEtc::default(); 1 << (shift + 1)], // kmer positions
            //kmp,
            contig: Vec::new(), // contig info
            repeat: AHashMap::new(),
            seed: seed,
        })
    }
    /*pub(crate) fn store_twobit(&mut self, pos: Position, b2: TwoBit) {
        self.b2[pos.byte_pos()] |= b2.pos_shift(pos).as_u8();
    }*/
    pub(crate) fn push_contig(&mut self, pos: Position, goffs: BasePos) {
        self.contig.push(Contig {
            twobit: pos,
            genomic: goffs,
        });
    }
    /// Adjust offset for contig.
    pub(crate) fn update_contig_genomic_offset(&mut self, goffs: BasePos) {
        if let Some(ctg) = self.contig.last_mut() {
            ctg.genomic = goffs;
        }
    }
    pub(crate) fn set_kmp(&mut self, mark: &XmerLoc) {
        dbg_print!("{} (stored)", mark);
        self.kmp[mark.idx].set(mark.p);
    }
    pub(crate) fn get_bitlen(&self) -> u8 {
        self.b2.len().trailing_zeros() as u8 + 2
    }
    pub(crate) fn is_on_last_contig(&self, pos: Position) -> bool {
        pos >= self.contig.last().unwrap().twobit
    }
    pub(crate) fn store_b2(&mut self, b2: TwoBit) -> Result<()> {
        if let Some(qb2) = self.b2.get_mut(self.b6pos >> 3) {
            *qb2 |= b2.as_u8() << (self.b6pos & 6);
            self.b6pos += 2;
            Ok(())
        } else {
            bail!("out of bounds for ks.b2: {}", self.get_pos())
        }
    }
    pub(crate) fn get_pos(&self) -> Position {
        Position::from_basepos((self.b6pos >> 1) as u64)
    }

    /// binary search contig lower boundary
    fn get_contig(&self, pos: Position) -> usize {
        let mut size = self.contig.len();
        let mut base = 0;
        while size > 0 {
            size /= 2;
            let mid = base + size;
            if self.contig[mid].twobit < pos {
                base = mid;
            }
        }
        base
    }
    pub(crate) fn get_contig_start_end_for_pos(&self, pos: Position) -> PosRange {
        let i = self.get_contig(pos);

        dbg_assert!(i < self.contig.len());
        dbg_assert!(self.contig[i].twobit <= pos);
        PosRange::from((self.get_twobit_before(i), self.get_twobit_after(i)))
    }
    fn get_twobit_after(&self, i: usize) -> Position {
        self.contig
            .get(i + 1)
            .map(|c| c.twobit)
            .unwrap_or(self.get_pos())
    }
    fn get_twobit_before(&self, i: usize) -> Position {
        i.checked_sub(1)
            .map(|t| self.contig[t].twobit)
            .unwrap_or(Position::default())
    }
    pub(crate) fn b2_for_pos(&self, pos: Position, is_repeat: bool) -> TwoBit {
        let byte_pos = pos.byte_pos();
        dbg_assert!(
            pos < self.get_pos() || is_repeat,
            "running into sequence head"
        );
        TwoBitx4::from(self.b2[byte_pos]).to_b2(pos, is_repeat)
    }
    pub(crate) fn extend_repetitive(&mut self, min_pos: Position, dist: BasePos) {
        let dist_u32 =
            u32::try_from(u64::from(dist)).expect("dist for repeat extension doesn't fit in u32");
        let repeat = self.repeat.entry(min_pos).or_insert((dist_u32, 0));
        repeat.1 = dist_u32;
    }

    pub(crate) fn replace_repetitive(
        &mut self,
        old_pos: Position,
        new_pos: Position,
        dist: BasePos,
    ) {
        let dist_u32 =
            u32::try_from(u64::from(dist)).expect("dist for repeat replacment doesn't fit in u32");
        if let Some(old_repeat) = self.repeat.remove(&old_pos) {
            let repeat = self.repeat.entry(new_pos).or_insert(old_repeat);
            repeat.1 += dist_u32;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rdbg::STAT_DB;
    #[test]
    fn it_works() -> Result<()> {
        //let mut ks = KmerStore::new(32); // allocates 16 gig of data
        let mut ks = KmerStore::new(2, 10_000, 0)?;
        let mut p = ExtPosEtc::default();
        let mut goffs = BasePos::default();

        ks.push_contig(p.pos(), BasePos::from(goffs));
        let mut stretch = BasePos::from(10_000);
        goffs += stretch;
        ks.update_contig_genomic_offset(goffs); // simulate N-stretch of 10000

        stretch = BasePos::from(64_u64);
        p = ExtPosEtc::from(p.pos() + Position::from(stretch)); // one readlength
        goffs += BasePos::from(p.pos());
        ks.push_contig(p.pos(), BasePos::from(goffs));
        ks.update_contig_genomic_offset(goffs);

        ks.push_contig(p.pos(), BasePos::from(goffs));
        stretch = BasePos::from(10_000);
        goffs += stretch;
        ks.update_contig_genomic_offset(goffs); // another N-stretch of 10000
        ks.push_contig(p.pos(), BasePos::from(goffs));

        let i = ks.get_contig(Position::from(BasePos::from(5_508_u64)));
        dbg_assert_eq!(ks.contig[i].twobit, Position::from(BasePos::from(64_u64)));
        Ok(())
    }
}
