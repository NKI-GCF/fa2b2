
#[macro_export]
macro_rules! dbgf {
	($l:literal) => ({
		if cfg!(debug_assertions) {
			eprint!("[{}:{}] {}\n", file!(), line!(), stringify!($l));
		}
		$l
	});
	($expr:expr, $fmt:literal$(, $opt:expr)*) => {
		match $expr {
			expr => {
				if cfg!(debug_assertions) {
					eprint!(concat!("[{}:{}] {} = ", $fmt, "\n"), file!(), line!(), stringify!($expr), &expr$(, $opt)*);
				}
				expr
			}
		}
	}
}
#[macro_export]
macro_rules! dbgx {
	($expr:expr) => ({
		if cfg!(debug_assertions) {
			dbg!($expr)
		} else {
			$expr
		}
	})
}

