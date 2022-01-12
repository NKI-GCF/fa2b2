//kmerstore
extern crate bincode;
extern crate flate2;
extern crate serde;

//extended_position
#[macro_use]
extern crate derive_new;
extern crate derive_more;

//kmer
extern crate num;
extern crate num_traits;
extern crate rand;

//for rdbg
#[macro_use]
extern crate lazy_static;

#[macro_use]
pub mod rdbg;
pub mod aln;
pub mod head_scope;
pub mod index;
pub mod kmerconst;
pub mod kmerstore;
pub mod mapping;
pub mod marker;
pub mod new_types;
pub mod past_scope;
pub mod scope;
