use crate::new_types::extended_position::ExtPosEtc;
use crate::rdbg::STAT_DB;
use std::{cmp, fmt};

#[derive(new, Copy, Clone, PartialEq, Eq)]
pub struct XmerLoc {
    pub idx: usize,
    pub p: ExtPosEtc,
}
impl XmerLoc {
    pub(crate) fn get(&self) -> Option<(usize, ExtPosEtc)> {
        if self.is_set() {
            Some((self.idx, self.p))
        } else {
            None
        }
    }
    pub(crate) fn get_idx(&self) -> usize {
        self.idx
    }
    pub(crate) fn is_set(&self) -> bool {
        self.idx != usize::MAX
    }
    /// kmers for unique positions for either strand within this mask are the scope of HeadScope,
    /// not the entire hash, but sufficiently large enough to account for repetitive sequences
    /// (transposons, simple or less simple repeats).
    pub(crate) fn get_scope_idx(&self, no_kmers: usize) -> usize {
        let mut scope_idx = self.p.as_basepos().as_usize() << 1;
        if !self.p.is_template() {
            scope_idx |= 1;
        }
        scope_idx & ((1 << no_kmers) - 1)
    }

    pub(crate) fn set(&mut self, idx: usize, p: ExtPosEtc) {
        // during rebuilding the strange case occurs that mark is not set, but p is (extension)
        let self_p_extension = self.p.extension();
        let p_extension = p.extension();
        dbg_assert!(self.is_set() || self.p.is_zero() || self_p_extension == p_extension);
        self.idx = idx;
        self.p = p;
    }
    pub(crate) fn get_thread_index(&self, bitlen: usize, ext_bits: usize) -> usize {
        self.idx >> (bitlen - ext_bits)
    }
}
impl Default for XmerLoc {
    fn default() -> Self {
        XmerLoc {
            idx: usize::MAX,
            p: Default::default(),
        }
    }
}

impl PartialOrd for XmerLoc {
    fn partial_cmp(&self, other: &XmerLoc) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for XmerLoc {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        match self.idx.cmp(&other.idx) {
            /*FIXME: make extenion count down, etc to simplify to this:
            cmp::Ordering::Equal => match self.p.cmp(&other.p) {
                cmp::Ordering::Equal => panic!(),
                x => x,
            },*/
            cmp::Ordering::Equal => match other.p.extension().cmp(&self.p.extension()) {
                cmp::Ordering::Equal => match self.p.pos().cmp(&other.p.pos()) {
                    cmp::Ordering::Equal => panic!(),
                    x => x,
                },
                x => x,
            },
            x => x,
        }
    }
}

impl fmt::Display for XmerLoc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{:x}] {}", self.idx, self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::new_types::position::BasePos;
    use rand::{thread_rng, Rng};

    impl XmerLoc {
        fn next(&mut self, ori: u64, is_template: bool) {
            if is_template {
                self.p.incr_pos()
            } else {
                self.p.decr_pos()
            }
            self.p.set_ori(ori & 1 != 0);
        }
    }
    #[test]
    fn forward() {
        let mut rng = thread_rng();
        let mut kl = XmerLoc::new(10, ExtPosEtc::from(BasePos::from(50_u64)));
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, true);
            dbg_assert_eq!(ori, if kl.p.is_template() { 0 } else { 1 });
        }
        dbg_assert_eq!(
            kl.p.as_u64() & !1,
            ExtPosEtc::from(BasePos::from(50_u64 + pick)).as_u64()
        );
    }
    #[test]
    fn reverse() {
        let mut rng = thread_rng();
        let mut kl = XmerLoc::new(10, ExtPosEtc::from(BasePos::from(50_u64)));
        let pick = rng.gen_range(20..50);
        for _ in 0..pick {
            let ori = rng.gen_range(0..2);
            kl.next(ori, false);
            dbg_assert_eq!(ori, if kl.p.is_template() { 0 } else { 1 });
        }
        dbg_assert_eq!(
            kl.p.as_u64() & !1,
            ExtPosEtc::from(BasePos::from(50_u64 - pick)).as_u64()
        );
    }
}
