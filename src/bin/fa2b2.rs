extern crate arrayvec;
extern crate bincode;
extern crate bio;
extern crate clap;
extern crate flate2;

extern crate fa2b2;

// target/release/fa2b2 -k 16 /net/NGSanalysis/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa test
//
// target/release/fa2b2 /home/roel/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz test

//use std::io::BufReader;
//use std::io::prelude::*;
//use flate2::bufread::MultiGzDecoder;
use clap::{App, Arg, SubCommand};

use anyhow::Result;
use fa2b2::aln::aln;
use fa2b2::index::index;

fn main() -> Result<()> {
    let matches = App::new("fa2b2")
        .version("1.0")
        .about("Read genome fasta file and write 2bit and region data ")
        .subcommand(
            SubCommand::with_name("index")
                .arg(
                    Arg::with_name("ref")
                        .short("r")
                        .long("ref")
                        .value_name("FASTA")
                        .help("The faidx'ed reference genome file, optionally bgzipped")
                        .index(1)
                        .required(true)
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("out")
                        .short("o")
                        .long("out")
                        .value_name("OUTPUT")
                        .help("The output file")
                        .index(2)
                        //.required(true)
                        .takes_value(true),
                )
                .arg(
                    Arg::with_name("repetition_max_dist")
                        .short("X")
                        .long("repetition-max-dist")
                        .value_name("repetition_max_dist")
                        .help("Maximum distance between kmers to be considered for repetition")
                        .takes_value(true),
                ),
        )
        .subcommand(
            SubCommand::with_name("aln").arg(
                Arg::with_name("ref")
                    .short("r")
                    .long("ref")
                    .value_name("FASTA")
                    .help("The faidx'ed reference genome file, optionally bgzipped")
                    .index(1)
                    .required(true)
                    .takes_value(true),
            ),
        )
        .get_matches();
    if let Some(matches) = matches.subcommand_matches("index") {
        index(matches)
    } else if let Some(matches) = matches.subcommand_matches("aln") {
        aln(matches)
    } else {
        unreachable!();
    }
}
