use crate::kmerconst::KmerConst;
use crate::kmerconst::XmerHash;
use crate::kmerstore::KmerStore;
use crate::new_types::{extended_position::EXT_MAX, position::Position};
use crate::rdbg::STAT_DB;
use crate::scope::Scope;
use crate::xmer_location::XmerLoc;
use anyhow::{anyhow, Result};
use bitvec::prelude::Lsb0;
use noodles_fastq as fastq;
use std::{cmp, fs, io};

pub struct Aligner<'a> {
    ks: KmerStore,
    kc: &'a KmerConst,
    scp: Scope<'a>,
}

impl<'a> Aligner<'a> {
    pub(crate) fn new(ks: KmerStore, kc: &'a KmerConst) -> Result<Self> {
        Ok(Aligner {
            ks,
            kc,
            scp: Scope::new(kc)?,
        })
    }
    // This function evaluates the number of multimappers.
    // An error indicates an xmer inconsistency: should be a mismatch to ref
    // in the window for this xmer. Read on until next xmer? else, try adjacent to picked median?
    fn se_mapping(&mut self, median_xmer: &XmerLoc) -> Result<()> {
        // binary search for dupbit 0 status
        let mut last_pos = Position::default();
        for x in 0..=EXT_MAX {
            let hash = self.kc.xmer_hash(median_xmer.idx, x);
            let test_p = self.ks.kmp[hash];
            match test_p.x().cmp(&x) {
                cmp::Ordering::Equal => {
                    let pos = test_p.pos();
                    // Because we haven't seen a non-DUPLICATE, we know there should follow more.
                    // x and hash lead to one baseindex. That should be position ordered.
                    dbg_assert!(last_pos < pos, "Xmer invalid by ref coordinate order.");

                    last_pos = pos;
                    if test_p.is_repetitive() {}
                    if test_p.is_last_on_ref() {
                        // no DUPLICATE: => last sequence that corresponded with this xmer.
                        break;
                    }
                }
                cmp::Ordering::Less => {
                    return Err(anyhow!("Xmer invalid by non-pos."));
                }
                cmp::Ordering::Greater => {}
            }
            // test_p.x() > x can happen: collisions
            // since the other extension was greater, this xmer & p got extended.
        }
        Ok(())
    }

    pub fn updated_optimal_xmers_only(&mut self, b: u8) -> Option<XmerLoc> {
        let b2 = match b {
            b'A' | b'a' => bitvec![u8, Lsb0; 0, 0],
            b'C' | b'c' => bitvec![u8, Lsb0; 1, 0],
            b'T' | b't' | b'U' | b'u' => bitvec![u8, Lsb0; 0, 1],
            b'G' | b'g' => bitvec![u8, Lsb0; 1, 1],
            n => {
                if n != b'N' && n != b'n' {
                    dbg_print!("Observed odd base {}, treating as ambiguous..", n);
                }
                // TODO / FIXME rather than excluding when Ns occur, try resolving those instances so
                // the seed selects against these as optima. Iterate over the possible sequences in place of the
                // ambiguous. Could try min..max. account which are the sweetspots. Weight is
                // the nr of positions this resolves.
                return None;
            }
        };
        //TODO: make sequence storage optional, to a separate file.
        if let Some(optimum_xmer) = self.scp.updated_median_xmer(&b2) {
            // FIXME: put the other strand in the idx top bits (saves reverse complementing).
            // Also, we need more alternative XmerLoc types for distinction along with traits !!
            return Some(optimum_xmer);
        }
        None
    }
    pub(crate) fn read_record(&mut self, record: &fastq::Record) -> Result<()> {
        for b in record.sequence().iter() {
            if let Some(median_xmer) = self.updated_optimal_xmers_only(*b) {
                self.se_mapping(&median_xmer)?;
            } else {
                //FIXME: N's ??
            }
        }
        Ok(())
    }
    pub(crate) fn align(
        &mut self,
        mut reader: fastq::Reader<io::BufReader<fs::File>>,
    ) -> Result<()> {
        let mut n = 0;
        for record in reader.records() {
            self.read_record(&record?)?;
            n += 1;
        }
        eprintln!("Processed {} reads", n);
        Ok(())
    }
}
