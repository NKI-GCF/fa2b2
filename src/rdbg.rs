use std::sync::Mutex;
use std::{
    cmp,
    collections::{HashMap, VecDeque},
};

lazy_static! {
    pub static ref STAT_DB: Mutex<StatDeq> = Mutex::new(StatDeq::new(500));
}

macro_rules! filelinestr {
    ($l:expr$(, $opt:expr)*) => ({
        concat!("[", file!(), ":", line!(), "] ", $l$(, $opt)*).to_string()
    });
}

/// formatted debug, included in dump but not in stats.
#[macro_export]
macro_rules! dbgf {
    ($l:literal) => ({
        if cfg!(debug_assertions) {
            STAT_DB.lock().unwrap().add(false, filelinestr!(stringify!($l)), stringify!($l));
            //eprintln!("[{}:{}] {}", file!(), line!(), stringify!($l));
        }
        $l
    });
    ($expr:expr, $fmt:literal$(, $opt:expr)*) => {
        match $expr {
            expr => {
                if cfg!(debug_assertions) {
                    STAT_DB.lock().unwrap().add(false, filelinestr!(stringify!($expr), " = ", $fmt),
                        format!(concat!("[{}:{}] {} => ", $fmt), file!(), line!(), stringify!($expr), &expr$(, $opt)*));
                    //eprintln!(concat!("[{}:{}] {} => ", $fmt), file!(), line!(), stringify!($expr), &expr$(, $opt)*);
                }
                expr
            }
        }
    }
}

/// begin a new round for logging (may trigger dump if dump_next() occurred)
#[macro_export]
macro_rules! dbg_restart {
    ($fmt:literal$(, $arg:expr)*) => ({
        if cfg!(debug_assertions) {
            STAT_DB.lock().unwrap().restart(filelinestr!(stringify!($fmt)), format!($fmt$(, $arg)*));
        }
    });
}

/// line for dump, not in stats
#[macro_export]
macro_rules! dbg_print {
    ($fmt:literal$(, $arg:expr)*) => ({
        if cfg!(debug_assertions) {
            STAT_DB.lock().unwrap().add(false, filelinestr!(stringify!($fmt)), format!($fmt$(, $arg)*));
        }
    });
}

/// line for dump, not in stats
#[macro_export]
macro_rules! dbg_print_if {
    ($cond:expr, $fmt:literal$(, $arg:expr)*) => ({
        if cfg!(debug_assertions) && $cond {
            STAT_DB.lock().unwrap().add(false, filelinestr!(stringify!($fmt)), format!($fmt$(, $arg)*));
        }
    });
}

/// for dump and in stats, conditionally dump_next()
#[macro_export]
macro_rules! dbg_dump_if {
    ($expr:expr, $cond:expr) => {{
        if let Some(mut db) = STAT_DB.lock().ok().filter(|_| cfg!(debug_assertions)) {
            match $expr {
                expr => {
                    db.add(
                        true,
                        format!(
                            "[{}:{}] {} => {}",
                            file!(),
                            line!(),
                            stringify!($expr),
                            &expr
                        ),
                        format!(
                            "[{}:{}] {} => {}",
                            file!(),
                            line!(),
                            stringify!($expr),
                            &expr
                        ),
                    );
                    if expr == $cond {
                        db.dump_next();
                    }
                    expr
                }
            }
        } else {
            $expr
        }
    }};
}

/// for dump and in stats
#[macro_export]
macro_rules! dbgx {
    ($expr:expr) => {{
        match $expr {
            expr => {
                if cfg!(debug_assertions) {
                    STAT_DB.lock().unwrap().add(
                        true,
                        format!(
                            "[{}:{}] {} => {:?}",
                            file!(),
                            line!(),
                            stringify!($expr),
                            &expr
                        ),
                        format!(
                            "[{}:{}] {} => {:?}",
                            file!(),
                            line!(),
                            stringify!($expr),
                            &expr
                        ),
                    );
                    //dbg!(&expr)
                }
                expr
            }
        }
    }};
}

///trigger dump directly.
#[macro_export]
macro_rules! dbg_dump {
    () => ({
        STAT_DB.lock().unwrap().dump(false);
        "--- end of dump ---\n"
    });
    ($fmt:literal$(, $arg:expr)*) => ({
        STAT_DB.lock().unwrap().dump(false);
        format!($fmt$(, $arg)*)
    });
}

#[macro_export]
macro_rules! dbg_assert {
    ($test:expr) => {{
        debug_assert!($test, "{}", dbg_dump!("Assert `{}' failed!", stringify!($test)));
    }};
    ($test:expr, $fmt:literal$(, $arg:expr)*) => {{
        debug_assert!($test, "{}", dbg_dump!($fmt$(, $arg)*));
    }};
}

#[macro_export]
macro_rules! dbg_assert_eq {
    ($a:expr, $b:expr) => {{
        debug_assert_eq!($a, $b, "{}", dbg_dump!("Assert `{}' == `{}' failed!", stringify!($a), stringify!($b)));
    }};
    ($a:expr, $b:expr, $fmt:literal$(, $arg:expr)*) => {{
        debug_assert_eq!($a, $b, "{}", dbg_dump!($fmt$(, $arg)*));
    }};
}

#[macro_export]
macro_rules! dbg_assert_ne {
    ($a:expr, $b:expr) => {{
        debug_assert_ne!($a, $b, "{}", dbg_dump!("Assert `{}' != `{}' failed!", stringify!($a), stringify!($b)));
    }};
    ($a:expr, $b:expr, $fmt:literal$(, $arg:expr)*) => {{
        debug_assert_ne!($a, $b, "{}", dbg_dump!($fmt$(, $arg)*));
    }};
}

#[macro_export]
macro_rules! dbg_panic {
    () => {{
        panic!("{}", dbg_dump!());
    }};
    ($fmt:literal$(, $arg:expr)*) => {{
        panic!("{}", dbg_dump!($fmt$(, $arg)*));
    }};
}

pub struct StatDeq {
    last: (String, String, bool),
    ct: usize,
    d: VecDeque<String>,
    h: HashMap<String, u32>,
    repeat: u32,
    dump_next: bool,
}

impl StatDeq {
    /// provide a message buffer of 'ct' lines, clamped between 100 and 10_000.
    pub fn new(ct: usize) -> Self {
        StatDeq {
            last: (String::new(), String::new(), false),
            ct: cmp::min(10000, cmp::max(ct, 100)),
            d: VecDeque::with_capacity(ct),
            h: HashMap::new(),
            repeat: 1,
            dump_next: false,
        }
    }
    pub fn restart(&mut self, fmt: String, msg: String) {
        if self.dump_next {
            self.dump(false);
            self.dump_next = false;
        }
        self.d.clear();
        self.add(false, fmt, msg);
    }
    pub fn add(&mut self, new_do_log: bool, fmt: String, msg: String) {
        let last_fmt = self.last.0.to_owned();
        let mut last_msg = self.last.1.to_owned();
        let do_log = self.last.2.to_owned();
        if last_fmt == fmt && last_msg == msg {
            self.repeat += 1;
        } else {
            if self.d.len() == self.ct {
                self.d.pop_front();
            }
            if self.repeat > 1 {
                last_msg.push_str(&format!(" ({}x)", self.repeat));
                self.repeat = 1;
            }
            self.d.push_back(last_msg);
            if do_log {
                *self.h.entry(last_fmt).or_insert(0) += 1;
            }
            self.last = (fmt, msg, new_do_log);
        }
    }
    pub fn dump_next(&mut self) {
        self.dump_next = true;
    }
    pub fn dump(&mut self, show_msg_counts: bool) {
        self.add(false, String::new(), String::new()); // to flush last message
        if show_msg_counts {
            eprintln!("--- Total file:line:message & counts: ---");
            let mut count_vec: Vec<_> = self.h.iter().collect();
            count_vec.sort_by(|a, b| a.0.cmp(b.0));
            for (msg, ct) in count_vec {
                eprintln!("{}\t{}", msg, ct);
            }
            eprintln!("--- Dump (last iteration): ---");
        }
        for msg in &self.d {
            eprintln!("{}", msg);
        }
    }
}
