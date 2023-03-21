extern crate arrayvec;
extern crate bincode;
extern crate clap;
extern crate flate2;

extern crate fa2b2;

// target/release/fa2b2 -k 16 /net/NGSanalysis/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa test
//
// target/release/fa2b2 /home/roel/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz test

//use std::io::BufReader;
//use std::io::prelude::*;
//use flate2::bufread::MultiGzDecoder;

use anyhow::Result;
use clap::{Parser, Subcommand};
use fa2b2::aln;
use fa2b2::index;
use std::path;

#[derive(Subcommand, Debug)]
enum Commands {
    /// indexes a reference file
    IndexCmd(index::IndexCmd),

    /// Align reads to a reference
    Aln(aln::AlnCmd),
}

/// Read genome fasta file and write 2bit and region data
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Fa2b2 {
    /// Sets a custom config file
    #[arg(short, long, value_name = "FILE")]
    config: Option<path::PathBuf>,

    /// Turn debugging information on
    #[arg(short, long, action = clap::ArgAction::Count)]
    debug: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

fn main() -> Result<()> {
    let fa2b2 = Fa2b2::parse();

    match fa2b2.command {
        Some(Commands::IndexCmd(index_cmd)) => index::index(index_cmd),
        Some(Commands::Aln(aln_cmd)) => aln::aln(aln_cmd),
        None => Ok(()),
    }
}
