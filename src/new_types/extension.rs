use crate::kmerloc::ExtPosEtc;

const EXT_SHIFT: u32 = 56;
const EXT_MASK: u64 = 0xFF00_0000_0000_0000;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Extension(u64);

impl Extension {
    pub fn as_u64(&self) -> u64 {
        self.0
    }
    pub fn max_value() -> Extension {
        Extension(EXT_MASK)
    }
}

// usize is assumed to be a base_extension, i.e. shifted to base.
impl From<usize> for Extension {
    fn from(base_ext: usize) -> Extension {
        let base_ext_u64 = u64::try_from(base_ext).unwrap();
        Extension(base_ext_u64.checked_shl(EXT_SHIFT).unwrap() & EXT_MASK)
    }
}

impl From<Extension> for usize {
    fn from(ext: Extension) -> usize {
        let base_ext_u64 = ext.0.checked_shr(EXT_SHIFT).unwrap();
        usize::try_from(base_ext_u64).unwrap()
    }
}

impl From<ExtPosEtc> for Extension {
    fn from(p: ExtPosEtc) -> Extension {
        Extension(p.as_u64() & EXT_MASK)
    }
}
