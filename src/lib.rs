//! Fa2b2 is a k-mer based DNA sequence indexer / mapper.
//!
//! Positions are stored in an array, k-mer based, but only for selected k-mers on the reference
//! sequence. Also, rather than using k-mers directly as index for position lookup in the array,
//! the k-mers are hashed to a certain extent.
//!
//! Due to increased hashing, the same k-mer can be stored for multiple positions. However, two
//! k-mers with distinct extension can also collide in the same hash. One has greater extension,
//! that one gets priority. The less extensively hashed k-mer gets extended further, resulting
//! in another different hash.
//!
//! Collisions are more work, but extension is pretty cheap. Actually a lot of time is spent in
//! selecting a median
//!
//! Differing k-mers cannot collide with the same extension in the same hash.
//!
//! Besides position, the hashing extent, and bits indicating strand, k-mer recurrence and the
//! presence of annotation.
//!
//!
//! The fa2b2 binary provides the two main subcommands index and aln.
//!

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
///
pub(crate) mod rdbg;
pub mod aln;
pub(crate) mod head_scope;
pub mod index;
pub(crate) mod kmerconst;
pub(crate) mod kmerstore;
pub(crate) mod mapping;
pub(crate) mod marker;
pub(crate) mod new_types;
pub(crate) mod past_scope;
pub(crate) mod scope;
