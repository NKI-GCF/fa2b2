//kmerstore
extern crate flate2;
extern crate rustc_serialize;
extern crate bincode;

//kmerconst
extern crate bit_reverse;

//kmerloc
#[macro_use]
extern crate derive_new;
//kmer
extern crate num;
extern crate num_traits;
extern crate rand;

//for rdbg
#[macro_use]
extern crate lazy_static;

#[macro_use]
pub mod rdbg;
pub mod kmer;
pub mod kmerstore;
pub mod kmerconst;
pub mod kmerloc;
pub mod occurrence;
pub mod marker;
