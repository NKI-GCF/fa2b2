
#[macro_export]
macro_rules! dbgf {
	($l:literal) => ({
		if cfg!(debug_assertions) {
			eprintln!("[{}:{}] {}", file!(), line!(), stringify!($l));
		}
		$l
	});
	($expr:expr, $fmt:literal$(, $opt:expr)*) => {
		match $expr {
			expr => {
				if cfg!(debug_assertions) {
					eprintln!(concat!("[{}:{}] {} = ", $fmt), file!(), line!(), stringify!($expr), &expr$(, $opt)*);
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

#[macro_export]
macro_rules! dbgln {
	($expr:expr) => ({
		if cfg!(debug_assertions) {
			eprintln!($expr);
		}
		$expr
	})
}

