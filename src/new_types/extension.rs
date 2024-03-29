// Roel Kluin, 2023, GPL v3

use crate::new_types::extended_position::{ExtPosEtc, EXT_MASK, EXT_SHIFT};
use crate::num::ToPrimitive;
use crate::rdbg::STAT_DB;

#[derive(Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Extension(u64);

impl Extension {
    pub(crate) fn as_u64(&self) -> u64 {
        self.0
    }
    pub(crate) fn max_value() -> Extension {
        Extension(EXT_MASK)
    }
}

// usize is assumed to be a base_extension, i.e. shifted to base.
impl From<usize> for Extension {
    fn from(base_ext: usize) -> Extension {
        let base_ext_u64 = base_ext.to_u64().unwrap();
        let masked = base_ext_u64 & (EXT_MASK >> EXT_SHIFT);
        dbg_assert_eq!(base_ext_u64, masked, "extension bleeds beyond!");
        Extension(base_ext_u64.checked_shl(EXT_SHIFT).unwrap() & EXT_MASK)
    }
}

impl From<Extension> for usize {
    fn from(ext: Extension) -> usize {
        ext.0.checked_shr(EXT_SHIFT).unwrap().to_usize().unwrap()
    }
}

impl From<ExtPosEtc> for Extension {
    fn from(p: ExtPosEtc) -> Extension {
        Extension(p.as_u64() & EXT_MASK)
    }
}
