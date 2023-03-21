use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::mapping::Aligner;
use anyhow::{ensure, Result};
use bincode::deserialize_from;
use clap::ArgMatches;
use noodles_fastq as fastq;
use std::{fs, io, path};

pub fn aln(matches: &ArgMatches) -> Result<()> {
    let fa_name = matches.value_of("ref").unwrap();
    let fq = matches.value_of("fastq1").unwrap();
    let reader = fs::File::open(fq)
        .map(io::BufReader::new)
        .map(fastq::Reader::new)?;

    let read_len = matches
        .value_of("read_length")
        .map(|v| v.parse())
        .transpose()?
        .unwrap();

    let ks_name = format!("{fa_name}.ks");
    ensure!(
        path::Path::new(&ks_name).exists(),
        "{ks_name} does not exist!"
    );
    let ks_file = io::BufReader::new(fs::File::create(ks_name)?);
    let ks: KmerStore = deserialize_from(ks_file)?;
    let bitlen = ks.get_bitlen();

    let kc = KmerConst::from_bitlen(bitlen, read_len, ks.seed);

    let mut aln = Aligner::new(ks, &kc)?;

    aln.align(reader)?;
    Ok(())
}
