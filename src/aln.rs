use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use anyhow::{ensure, Result};
use bincode::deserialize_from;
use clap::ArgMatches;
use noodles_fastq as fastq;
use std::{fs, io, path};

pub fn aln(matches: &ArgMatches) -> Result<()> {
    let fa_name = matches.value_of("ref").unwrap();
    let fq = matches.value_of("fastq1").unwrap();
    let mut reader = fs::File::open(fq)
        .map(io::BufReader::new)
        .map(fastq::Reader::new)?;

    let ks_name = format!("{}.ks", fa_name);
    ensure!(
        path::Path::new(&ks_name).exists(),
        "{} does not exist!",
        ks_name
    );
    let ks_file = io::BufReader::new(fs::File::create(ks_name)?);
    let ks: KmerStore<u64> = deserialize_from(ks_file)?;
    let bitlen: usize = ks.get_bitlen();

    let kc = KmerConst::from_bitlen(bitlen);

    let mut n = 0;

    for res in reader.records() {
        let record = res?;
        n += 1;
    }

    Ok(())
}
