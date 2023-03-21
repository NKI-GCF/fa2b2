use crate::kmerconst::KmerConst;
use crate::kmerstore::KmerStore;
use crate::mapping::Aligner;
use anyhow::{ensure, Result};
use bincode::deserialize_from;
use clap::Args;
use noodles_fastq as fastq;
use std::{fs, io, path};

#[derive(Args, Debug)]
pub struct AlnCmd {
    /// The faidx'ed reference genome file, optionally bgzipped
    #[arg(short, long, value_name = "FASTA", required = true)]
    ref_file: path::PathBuf,

    /// The forward fastq file
    #[arg(short, long, value_name = "FASTQ", required = true)]
    fq1: path::PathBuf,

    /// The reverse fastq file if sequeced paired-end
    #[arg(short, long, value_name = "FASTQ")]
    opt_fq2: Option<path::PathBuf>,

    /// Length of sequence reads
    #[arg(short, long, value_name = "read_length", required = true)]
    read_len: u16,
}

pub fn aln(cmd: AlnCmd) -> Result<()> {
    let reader = fs::File::open(cmd.fq1)
        .map(io::BufReader::new)
        .map(fastq::Reader::new)?;

    let ks_name = format!("{:?}.ks", cmd.ref_file);
    ensure!(
        path::Path::new(&ks_name).exists(),
        "{ks_name} does not exist!"
    );
    let ks_file = io::BufReader::new(fs::File::create(ks_name)?);
    let ks: KmerStore = deserialize_from(ks_file)?;
    let bitlen = ks.get_bitlen();

    let kc = KmerConst::from_bitlen(bitlen, cmd.read_len, ks.seed);

    let mut aln = Aligner::new(ks, &kc)?;

    aln.align(reader)?;
    Ok(())
}
