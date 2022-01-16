use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::new_types::extended_position::ExtPosEtc;
use crate::new_types::{
    position::{BasePos, Position},
    twobit::TwoBit,
    xmer::Xmer,
};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::xmer_location::XmerLoc;
use anyhow::Result;
use smallvec::{smallvec, SmallVec};
use std::fmt;
use to_default::ToDefault;
const SCOPE_WIDTH_MAX: usize = 128;

pub struct HeadScope<'a> {
    kc: &'a KmerConst,
    p: ExtPosEtc,
    d: SmallVec<[Xmer; SCOPE_WIDTH_MAX]>,
    z: SmallVec<[usize; SCOPE_WIDTH_MAX]>,
    mark: XmerLoc,
    i: usize, //TODO: increment this only per mod_i, then could be u32, could matter for padding
    mod_i: usize,
    pub(crate) repetitive: u32,
    n_stretch: BasePos,
    goffs: BasePos,
    period: Position,
}

impl<'a> HeadScope<'a> {
    pub(crate) fn new(kc: &'a KmerConst) -> Self {
        HeadScope {
            kc,
            p: ExtPosEtc::default(),
            d: smallvec![Xmer::new(); kc.no_kmers],
            z: (0..kc.no_kmers).into_iter().collect::<SmallVec<_>>(),
            mark: XmerLoc::new(usize::max_value(), ExtPosEtc::default()),
            i: 0,
            mod_i: 0,
            repetitive: 0,
            n_stretch: BasePos::default(),
            goffs: BasePos::default(),
            period: Position::default(),
        }
    }

    // .i & .p increments en kmer .d[] updates vinden plaats.
    pub(crate) fn complete_and_update_mark(
        &mut self,
        b2: TwoBit,
        ks: &mut KmerStore,
    ) -> Result<()> {
        self.manage_former_ns_and_period(ks, b2)?;
        if self.update() {
            // one mark is added, and one leaving. both influence mark (and order).
            let i = self.pick_mark();
            if self.d[i].pos == self.mark.p.pos() {
                return Ok(());
            }
            let mut mark = self.d[i].get_hash_and_p(self.kc, 0);
            self.set_mark(&mark);
            let orig_pos = mark.p.pos();
            // this seems to be hotlooping
            while self.try_store_mark(ks, &mut mark)? {
                if mark.p.pos() != orig_pos {
                    if self.kc.extend_xmer(&mut mark).is_ok() {
                        // extending some pase baseidx. TODO: if frequently the same recurs,
                        // it might be worthwhile to store the reverse complement in a temp
                        // we should not store a past index in self.mark.p !!
                    } else {
                        dbg_print!("couldn't extend: {} ..?", mark);
                        break;
                    }
                } else {
                    if mark.p.extend().is_err() {
                        break;
                    }
                    mark = self.d[i].get_hash_and_p(self.kc, mark.p.x());
                    self.set_mark(&mark);
                }
            }
        }
        self.increment(b2);
        Ok(())
    }

    fn manage_former_ns_and_period(&mut self, ks: &mut KmerStore, b2: TwoBit) -> Result<()> {
        self.finalize_n_stretch(ks);
        if self.period.is_set() {
            let pd = self.period;
            let pos = self.p.pos();
            dbg_assert!(pd <= pos, "{} {}", pd, self.p);
            if ks.b2_for_p(pos - pd, true)? == b2 {
                self.update_repetitive(ks, pd);
            } else {
                self.period = Position::default();
            }
        }
        Ok(())
    }
    /// Finalize an N-stretch, update stored offsets and prepare to process sequence.
    pub(crate) fn finalize_n_stretch(&mut self, ks: &mut KmerStore) {
        if self.n_stretch.is_set() {
            dbg_print!("added new contig. Ns:{}", self.n_stretch);
            self.goffs += ks.offset_contig(&mut self.n_stretch);

            // clear all to rebuild at the start of a new contig.
            self.p.clear_extension();
            self.mod_i = 0;
            // If a repetition ends in an N-stretch, thereafter offset to period
            // may differ or the repetition could be different or entirely gone.
            // TODO: allow repetition to include N-stretch - if both sides of N-stretch show the same repetition.
            self.period = Position::default();
        }
    }

    pub(crate) fn reset_for_new_contig(&mut self) -> Position {
        // we start with no offset on contig, if starting with N's, the stored goffs gets updated
        self.goffs = BasePos::default();
        self.n_stretch.to_default();
        self.repetitive = 0;
        self.p.pos()
    }
    pub(crate) fn elongate_n_stretch(&mut self, ks: &mut KmerStore, pos: Position) {
        if self.i != 0 {
            dbg_print!("started N-stretch at {}.", self.p.basepos());
            self.goffs.add_assign(self.i as u64);
            ks.push_contig(pos, self.goffs);

            self.n_stretch.to_default();
            self.i = 0;
        }
        self.n_stretch.add_assign(1_u64);
    }

    fn update_repetitive(&mut self, ks: &mut KmerStore, pd: Position) {
        // XXX waarom werkt dit op mark??
        self.repetitive += 1;
        let idx = self.mark.get_idx();
        let stored = ks.kmp[idx];
        if stored.is_set() {
            let mark_pos = self.mark.p.pos();
            let stored_pos = stored.pos();
            if mark_pos != stored_pos {
                if let Some(dist) = mark_pos.get_if_mark_on_period(stored_pos, pd) {
                    if mark_pos > stored_pos {
                        ks.extend_repetitive(mark_pos, dist);
                    } else {
                        dbg_print!(
                            "repetitive occurs before already stored [{:x}] p {} <=> stored {}",
                            idx,
                            self.mark.p,
                            stored
                        );
                        let _x = self.p.x();
                        dbg_assert!(_x > 0);
                        ks.replace_repetitive(stored_pos, mark_pos, dist);
                        // FIXME: moet positie nu niet gezet worden, bij less?
                    }
                }
            }
        } else {
            // XXX this occurs, is it an edge case or a bug?
            dbg_print!("repeat with unset mark.idx (b2 corresponds) ??");
        }
    }

    fn try_store_mark(&mut self, ks: &mut KmerStore, mark: &mut XmerLoc) -> Result<bool> {
        if self.period.is_set() && ks.kmp[mark.idx].is_set() {
            ks.kmp[mark.idx].set_repetitive();
            return Ok(false);
        }
        let old_stored_p = ks.kmp[mark.idx];

        if old_stored_p.is_zero() {
            ks.set_kmp(&mark);
            return Ok(false);
        }
        if old_stored_p.is_replaceable_by(mark.p) {
            if old_stored_p.pos() == mark.p.pos() {
                // set and already mark.p. Leave the bit states.
                return Ok(false);
            }
            ks.set_kmp(&mark);
            if old_stored_p.x() == mark.p.x() {
                // same extension means same base k-mer origin. this is a duplicate.
                ks.kmp[mark.idx].mark_more_recurs_upseq();
            }
            dbg_print!("{} -> ?", mark);
            mark.p = old_stored_p; // to be extended next
        } else if old_stored_p.extension() == mark.p.extension() {
            // Note: same extension and hash means same k-mer origin: identical k-mer sequence.

            // If a kmer occurs multiple times within an extending readlength (repetition),
            // only the first gets a position. During mapping this should be kept in mind.
            if let Some(dist) = self.dist_if_repetitive(ks, old_stored_p, mark.p) {
                dbg_assert!(dist < self.p.pos() || self.p.is_zero());
                self.period = dist;
                ks.kmp[mark.idx].set_repetitive();
                return Ok(false);
            }
            ks.kmp[mark.idx].mark_more_recurs_upseq();
        }
        // collision between hashes, the one in mark.p will be extended and tried again.
        Ok(true)
    }
}

impl<'a> Scope for HeadScope<'a> {
    fn get_pos(&self) -> Position {
        self.p.pos()
    }

    /// add twobit to k-mers, update k-mer vec, increment pos and update orientation
    /// true if we have at least one kmer.
    fn update(&mut self) -> bool {
        // XXX: function is hot
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.mod_i];
            self.mod_i += 1;
            if self.mod_i == self.kc.no_kmers {
                self.mod_i = 0;
            }
            self.d[self.mod_i] = old_d;
            self.d[self.mod_i].pos = self.p.pos();
            true
        } else {
            false
        }
    }
    fn pick_mark(&mut self) -> usize {
        let med = self.kc.no_kmers >> 1;
        let i = self
            .z
            .select_nth_unstable_by(med, |&a, &b| self.d[a].cmp(&self.d[b]))
            .1;
        *i
    }

    fn dist_if_repetitive(
        &self,
        ks: &KmerStore,
        stored_p: ExtPosEtc,
        min_p: ExtPosEtc,
    ) -> Option<Position> {
        // FIXME: the contig check was removed. for the upper bound that makes sense for head
        // scope, as upperbound is pos_max rather than previous
        let stored_pos = stored_p.pos();
        let mark_pos = min_p.pos();
        dbg_assert!(mark_pos > stored_pos);
        if ks.is_on_last_contig(stored_pos) {
            let dist = mark_pos - stored_pos;
            if dist < ks.rep_max_dist {
                return Some(dist);
            }
        }
        None
    }

    fn set_mark(&mut self, mark: &XmerLoc) {
        dbg_print!("{} (mark, head)", mark);
        self.mark = *mark;
    }

    fn increment(&mut self, b2: TwoBit) {
        // first bit is strand bit, set according to kmer orientation bit.
        self.p.set_ori(self.d[self.mod_i].update(self.kc, b2));
        self.p.incr_pos();
        self.i += 1;
    }
}

impl<'a> fmt::Display for HeadScope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.unshift_pos() as usize;
        let mp = self.mark.p.unshift_pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat((p - self.kc.read_len) * 5);
        let r = p - mp;
        if r == 0 {
            let x = self.kc.read_len - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.get_idx(), n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.read_len {
            let x = self.kc.read_len - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.get_idx(), n << 2, o, s)
        } else {
            //let l = self.kc.read_len - r - n;
            //let ls = if o {" ".repeat(o) + "|"} else {String::from("")};
            //let rs = if l {String::from("|") + &" ".repeat(l << 2)} else {String::from("")};
            //write!(f, "{2}<{3}|{: ^1$x}|{4}>", self.mark.idx, n << 2, o, ls, rs)
            write!(f, "")
        }
    }
}
