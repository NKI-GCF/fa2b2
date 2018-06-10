extern crate bio;
extern crate flate2;
extern crate arrayvec;
extern crate rustc_serialize;
extern crate bincode;
extern crate clap;

// target/release/fa2b2 -k 16 /net/NGSanalysis/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa test
//
// target/release/fa2b2 -k 16 /home/roel/Downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz test

use std::io::{self, Write};
use bio::io::fasta::IndexedReader;

use clap::{App,Arg};

mod kmer;
use kmer::marker::KmerMarker;

fn main() {
    let matches = App::new("fa2b2")
                .version("1.0")
                .about("Read genome fasta file and write 2bit and region data ")
                .arg(Arg::with_name("ref")
                        .short("r")
                        .long("ref")
                        .value_name("FASTA")
                        .help("The .fa reference genome file")
                        .index(1)
                        .required(true)
                        .takes_value(true))
                .arg(Arg::with_name("out")
                        .short("o")
                        .long("out")
                        .value_name("OUTPUT")
                        .help("The output file")
                        .index(2)
                        .required(true)
                        .takes_value(true))
                .get_matches();

        let fasta = matches.value_of("ref").unwrap();

        let mut idxr = IndexedReader::from_file(&fasta).unwrap_or_else(|_| panic!("Error opening reference genome"));
        let chrs = idxr.index.sequences();

        let outfile = matches.value_of("out").unwrap();

        let genomesize: usize = chrs.iter().map(|x| x.len as usize).sum();

        let mut kmermarker = KmerMarker::new(genomesize);

        let mut p = 0;
        println!("Chromosome\trunning unique count");
        for chr in &chrs {
            let mut seq = Vec::with_capacity(chr.len as usize);
            idxr.fetch_all(&chr.name).expect(&format!("Error fetching {}.", &chr.name));
            idxr.read(&mut seq).expect(&format!("Error reading {}.", &chr.name));

            p += kmermarker.markcontig(&seq, p);
            println!("{}\t{}", &chr.name, p);
            io::stdout().flush().unwrap();
        }

        println!("Writing first occurances per kmer to disk");
        kmermarker.write(outfile).unwrap();
}


