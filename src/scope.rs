use crate::kmer::Kmer;
use crate::kmerconst::KmerConst;
use crate::kmerloc::{KmerLoc, PriExtPosOri};
use crate::kmerstore::KmerStore;
use crate::rdbg::STAT_DB;
use anyhow::{ensure, Result};
use std::{cmp, fmt};

pub struct Scope<'a> {
    pub kc: &'a KmerConst,
    pub p: u64,
    d: Vec<Kmer<u64>>, // misschien is deze on the fly uit ks te bepalen?
    pub mark: KmerLoc<u64>,
    pub i: usize,
    pub mod_i: usize,
    pub plim: (u64, u64),
    pub period: u64,
}

impl<'a> Scope<'a> {
    pub fn new(plim: (u64, u64), kc: &'a KmerConst, p: u64) -> Self {
        Scope {
            kc,
            p,
            d: vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers],
            mark: KmerLoc::new(usize::max_value(), p.extension()),
            i: 0,
            mod_i: 0,
            plim,
            period: 0,
        }
    }

    pub fn rebuild<T: PriExtPosOri + fmt::LowerHex>(
        ks: &KmerStore<T>,
        kc: &'a KmerConst,
        p: &T,
        idx: usize,
    ) -> Result<Self> {
        /* Er gaat iets mis, blijkt uit repetitie telling die afwijkt afhankelijk
        van de scope heropbouw. Indien met plim.0, geeft dit andere repetition
        aantallen dan wanneer de eerste 2bit in scope wordt gebruikt. Voor de
        meeste contigs meer repetitie, enkele minder. Duizendtallen op de primary
        contigs.
        De positie waarvoor herbouwt wordt zou een repetition of duplicaat kunnen
        zijn, maar repetition wordt alleen geteld bij nieuwe sequentie. Dus afh
        van scope rebuilding moeten er verschillende xmers 'bezet' zijn, wat dan
        invloed heeft op de telling bij nieuwe sequence, lijkt mij de logischte
        verklaring. */

        // FIXME: dit zou hetzelfde resultaat moeten geven!
        // TODO: verbeter rebuild code en test impact hierop.
        // Minder of geen verschil bij deze twee alternatieven is beter:
        // self.p = p.extension() | plim.0;
        let pos = p.pos();
        ensure!(pos != PriExtPosOri::no_pos());
        let plim = ks.get_contig_start_end_for_p(pos);
        let bound = kc.get_kmer_boundaries(pos, plim);

        let mut scp = Scope {
            kc,
            p: bound.0 | p.extension(),
            d: vec![Kmer::new(kc.kmerlen as u32); kc.no_kmers],
            mark: KmerLoc::new(usize::max_value(), p.extension()),
            i: 0,
            mod_i: 0,
            plim,
            period: 0,
        };

        loop {
            let b2 = ks.b2_for_p(scp.p).unwrap();
            dbg_print!("=> b2 {:x}, p: {:#x} <=", b2, scp.p);
            if scp.increment(b2) {
                // we weten extension op voorhand.
                let oks = Some(ks);
                if (scp.set_if_optimum(0, scp.mark.p.x(), oks) && scp.all_kmers())
                    || scp.mark_set_considering_leaving(oks)?
                {
                    if p.same_pos_and_ext(scp.mark.p) {
                        break;
                    }
                    if scp.mark.get_idx() == idx {
                        dbg_print!(
                            "idx {:x} observed but for {:#x}, not {:#x}",
                            idx,
                            scp.mark.p,
                            p
                        );
                    }
                    // XXX ik zou een assertion hier logischer vinden
                    /*assert!(
                        scp.p.pos() < bound.1,
                        "kmer {:x} not observed for {:x} !!",
                        idx,
                        p
                    );*/
                    if scp.p.pos() >= bound.1 {
                        dbg_print!("kmer {:x} not observed for {:x} !!", idx, p);
                        scp.p.clear();
                        break;
                    }
                }
            }
        }
        return Ok(scp);
    }

    /// return toont of er een oplossing was
    pub fn handle_mark<T>(&mut self, ks: &mut KmerStore<T>) -> Result<bool>
    where
        T: PriExtPosOri + fmt::LowerHex + Copy,
    {
        let (min_idx, min_p) = self.mark.get();
        let stored_p = ks.kmp[min_idx];

        if stored_p.is_replaceable_by(min_p) {
            //dbg_print!("[{:#x}] (={:#x}) <= {:#x}", min_idx, stored_p, min_p);
            // dbg_print!("{}", self);
            if dbgx!(stored_p.is_set_and_not(min_p)) {
                let mut new_scp = Scope::rebuild(ks, self.kc, &stored_p, min_idx)?;
                if !new_scp.p.is_set() {
                    // unset when kmer is not observed (zie dbg_print..)
                    return Ok(false);
                }
                if let Some(storage) = ks.kmp.get_mut(min_idx) {
                    storage.set(min_p);
                }
                dbg_print!("resolving past for [{:x}], {:#x}", min_idx, stored_p);
                if new_scp.extension_gave_mark(ks)? && !new_scp.handle_mark(ks)? {
                    dbg_print!("unresolved new_scp mark");
                }
            } else if ks.kmp[min_idx].is_no_pos() {
                ks.kmp[min_idx].set(min_p);
                // If already set it was min_p. Then leave dupbit state.
            }
            return Ok(true);
        }
        if stored_p.extension() == min_p.extension() {
            // If a kmer occurs multiple times within an extending readlength (repetition),
            // only the first gets a position. During mapping this rule should also apply.
            // Mark / store recurring xmers. Skippable if repetitive on contig. Else mark as dup.
            let stored_pos = stored_p.pos();
            assert!(self.mark.p.pos() > stored_pos);
            if dbgx!(stored_pos >= self.plim.0.pos() && stored_pos < self.plim.1.pos()) {
                let dist = self.mark.p.pos() - stored_pos;
                if dbgx!(dist < self.kc.repetition_max_dist) {
                    self.period = dist;
                    // niet gezet, doen we niet bij repetition, extensie is niet nodig.
                    return Ok(false);
                }
            }
            ks.kmp[min_idx].set_dup();
        } else {
            dbg_print!("\t\t<!>");
        }
        // in repetitive dna, mark may not be assigned
        self.extension_gave_mark(ks)
    }

    fn can_extend(&self) -> bool {
        self.p.x() + 1 < self.kc.extent.len()
    }

    // reads b2 for extension, but we cannot add, e.g. running into end of contig => false.
    fn extend<T: PriExtPosOri>(&mut self, ks: &KmerStore<T>) -> Result<()> {
        self.p.extend();
        dbg_assert!(self.p.pos() >= self.mark.p.pos());

        if self.p.pos() - self.mark.p.pos() > (self.kc.afstand(self.p.x()) * 2) as u64 {
            let b2 = ks.b2_for_p(self.p).unwrap();
            let has_kmer = self.increment(b2);
            assert!(has_kmer);
            ensure!(!self.is_p_beyond_contig(), "unresolved past extension");
        }
        Ok(())
    }

    /// true: too few twobits or at contig end; true: extended and mark set.
    fn extension_gave_mark<T: PriExtPosOri>(&mut self, ks: &KmerStore<T>) -> Result<bool> {
        while self.can_extend() {
            if let Err(e) = self.extend(ks) {
                dbg_print!("{}", e);
                self.p.clear();
                break;
            }
            if self.set_next_mark(Some(ks))? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn is_xmer_complete(&self, x: usize) -> bool {
        self.i >= self.kc.kmerlen + self.kc.afstand(x)
    }

    fn is_p_beyond_contig(&self) -> bool {
        dbgx!(self.p.pos() >= self.plim.1)
    }

    /// add twobit to k-mers, update k-mer vec, increment pos and update orientation
    /// true if we have one kmer.
    pub fn increment(&mut self, b2: u8) -> bool {
        // XXX: function is hot
        if self.i >= self.kc.kmerlen {
            let old_d = self.d[self.mod_i];
            self.mod_i += 1;
            if self.mod_i == self.kc.no_kmers {
                self.mod_i = 0;
            }
            self.d[self.mod_i] = old_d;
        }
        // first bit is strand bit, set according to kmer orientation bit.
        self.p &= !1;
        self.p += 2 + self.d[self.mod_i].update(b2);
        self.i += 1;
        self.i >= self.kc.kmerlen
    }

    /// sufficient kmers to have mimimum and do we have minimum?
    fn mark_set_considering_leaving<T: PriExtPosOri>(
        &mut self,
        oks: Option<&KmerStore<T>>,
    ) -> Result<bool> {
        if self.all_kmers() && self.mark.is_set() {
            if self.is_mark_leaving() {
                return self.set_next_mark(oks);
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }

    // hierin vinden .i & .p increments and kmer .d[] update plaats.
    // returns true if all kmers observed and one is set as mark.
    pub fn complete_and_update_mark<T: PriExtPosOri>(
        &mut self,
        b2: u8,
        oks: Option<&KmerStore<T>>,
    ) -> Result<bool> {
        // XXX: function is very hot
        if self.increment(b2) {
            // mark.p hoeft niet gezet te zijn.
            for x in 0..=self.p.x() {
                if self.set_if_optimum(0, x, oks) {
                    // als we een nieuwe xmer tegenkwamen, wel genoeg kmers?
                    return Ok(self.all_kmers());
                }
            }
            self.mark_set_considering_leaving(oks)
        } else {
            Ok(false)
        }
    }

    /// is occ complete? call na complete_kmer() - self.i increment.
    fn all_kmers(&self) -> bool {
        self.i >= self.kc.venster
    }

    fn is_mark_leaving(&self) -> bool {
        self.p.pos() >= self.mark.p.pos() + (self.kc.no_xmers(self.p.x()) << 1) as u64
    }

    /// bepaal van alle kmers het nieuwe minimum/optimum (na leaving mark of extensie)
    /// geen enkele optimum kan voorkoment bij repetitive DNA
    fn set_next_mark<T: PriExtPosOri>(&mut self, oks: Option<&KmerStore<T>>) -> Result<bool> {
        self.mark.reset();
        let x = self.p.x();
        // reverse is logischer omdat we bij gelijke extensie voor een lagere positie kiezen.
        for i in (0..self.kc.no_xmers(x)).rev() {
            let _ = self.set_if_optimum(dbgx!(i), x, oks);
        }
        Ok(self.mark.is_set())
    }

    /// voor een offset i en extensie x, maak de kmer/hash en zet mark + return true als optimum.
    fn set_if_optimum<T: PriExtPosOri>(
        &mut self,
        i: usize,
        x: usize,
        oks: Option<&KmerStore<T>>,
    ) -> bool {
        // XXX function is hot
        if self.is_xmer_complete(x) {
            let e = self.kc.get_kmers(x);
            let base = self.i - self.kc.kmerlen;
            let d = i + e.0 as usize;
            let nk = self.kc.no_kmers;
            let kmer1 = self.d[base.wrapping_sub(d) % nk];
            let mut hash = kmer1.get_idx(e.0 <= e.1);
            let mut p = match hash.cmp(&self.mark.get_idx()) {
                cmp::Ordering::Less => self.p.with_ext(x) - (d << 1) as u64,
                cmp::Ordering::Greater => return false,
                cmp::Ordering::Equal => {
                    let p = self.p.with_ext(x) - (d << 1) as u64;
                    if p >= self.mark.p {
                        return false;
                    }
                    p
                }
            };
            p ^= match e.0.cmp(&e.1) {
                cmp::Ordering::Less => {
                    hash ^= self.d[base.wrapping_sub(i + e.1 as usize) % nk].get_idx(true);
                    p & 1
                }
                cmp::Ordering::Greater => {
                    hash ^= self.d[base.wrapping_sub(i + e.1 as usize) % nk].get_idx(false);
                    !p & 1
                }
                cmp::Ordering::Equal => (p ^ kmer1.dna) & 1,
            };
            if let Some(rep) = oks.and_then(|ks| ks.repeat.get(&hash)) {
                let ks = oks.unwrap();
                let stored_p = &ks.kmp[hash];
                if stored_p.is_set()
                    && p.pos() != stored_p.pos()
                    && dbgf!(
                        p.pos() < stored_p.pos() + rep.1 as u64,
                        "{} {:#x} Rep:{:?}",
                        hash,
                        rep
                    )
                {
                    return false;
                }
            }
            dbgf!(
                self.mark.set(hash, p, x),
                "{:?} [{:x}] = {:x} | x({})",
                hash,
                p,
                x
            );
            // XXX maybe return self.all_kmers() ??
            true
        } else {
            false
        }
    }
}

impl<'a> fmt::Display for Scope<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = self.p.pos() as usize;
        let mp = self.mark.p.pos() as usize;
        let n = self.kc.kmerlen + self.p.x();
        let o = " ".repeat(((p >> 1) - self.kc.venster) * 5);
        let r = (p - mp) >> 1;
        if r == 0 {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                " ".repeat(x << 2) + "|"
            } else {
                String::from("")
            };
            write!(f, "{2}<{3}{: ^1$x}>", self.mark.get_idx(), n << 2, o, s)
        } else if r + self.kc.kmerlen == self.kc.venster {
            let x = self.kc.venster - n;
            let s = if x != 0 {
                String::from("|") + &" ".repeat(x << 2)
            } else {
                String::from("")
            };
            write!(f, "{2}<{: ^1$x}{3}>", self.mark.get_idx(), n << 2, o, s)
        } else {
            //let l = self.kc.venster - r - n;
            //let ls = if o {" ".repeat(o) + "|"} else {String::from("")};
            //let rs = if l {String::from("|") + &" ".repeat(l << 2)} else {String::from("")};
            //write!(f, "{2}<{3}|{: ^1$x}|{4}>", self.mark.idx, n << 2, o, ls, rs)
            write!(f, "")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl<'a> Scope<'a> {
        /// de current kmer.
        fn kmer(&self) -> Kmer<u64> {
            self.d[self.mod_i]
        }

        /// add b2 to kmer, move .p according to occ orientation
        /// return whether required nr of kmers for struct scope were seen.
        fn update_mark<T: PriExtPosOri>(&mut self, x_start: usize, oks: Option<&KmerStore<T>>) {
            if self.i >= self.kc.kmerlen {
                for x in x_start..=self.p.x() {
                    if self.set_if_optimum(0, x, oks) {
                        return;
                    }
                }
            }
        }
    }

    const SEQLEN: usize = 250;

    #[test]
    fn test_push_b2() {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut occ = Scope::new((0, 100), &kc, 0);
        for _ in 0..6 {
            let _ = occ.increment(1);
            occ.update_mark::<u64>(0, None);
        }
        let mut kmer = occ.kmer();
        dbg_assert_eq!(kmer.dna, 0x55);
        dbg_assert_eq!(kmer.rc, 0xff);
        dbg_assert_eq!(occ.p, 0xd);

        occ.increment(2);
        occ.update_mark::<u64>(0, None);
        kmer = occ.kmer();
        dbg_assert_eq!(kmer.dna, 0x95);
        dbg_assert_eq!(kmer.rc, 0xfc);
        dbg_assert_eq!(occ.p, 0xf);
        dbg_assert_eq!(kmer.get_idx(true), 0x6a);
    }
    #[test]
    fn occurrence() {
        let kc = KmerConst::new(SEQLEN, 1000);
        let mut occ = Scope::new((0, 100), &kc, 0);
        for _ in 0..8 {
            let _ = occ.increment(0);
            occ.update_mark::<u64>(0, None);
        }
        for _ in 0..8 {
            let _ = occ.increment(1);
            occ.update_mark::<u64>(0, None);
        }
        for _ in 0..8 {
            let _ = occ.increment(2);
            occ.update_mark::<u64>(0, None);
        }
        for _ in 0..8 {
            let _ = occ.increment(3);
            occ.update_mark::<u64>(0, None);
        }
    }
}
