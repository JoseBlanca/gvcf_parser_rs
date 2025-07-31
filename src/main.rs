use clap::Parser;
use polars::prelude::ParquetWriter;
use polars::prelude::*;
use std::fs::File;
use std::path::Path;
use std::path::PathBuf;

use gvcfparser::errors::VcfParseError;
use gvcfparser::gvcf_parser::{GVcfRecord, GVcfRecordIterator, VcfResult};

/// Extract variant regions from a gVCF and save them to a Parquet file.
#[derive(Parser, Debug)]
#[command(name = "gvcf_to_parquet")]
#[command(author = "Jose Blanca")]
#[command(version = "0.1.0")]
#[command(about = "Extracts variant spans from a .g.vcf.gz and stores them in a Parquet file.", long_about = None)]
struct Args {
    /// Input .g.vcf.gz path
    #[arg(short, long)]
    input: PathBuf,

    /// Output Parquet file path
    #[arg(short, long)]
    output: PathBuf,
}

pub fn save_var_regions_as_parquet<I, P>(iterator: I, output_path: P) -> PolarsResult<()>
where
    I: Iterator<Item = VcfResult<GVcfRecord>>,
    P: AsRef<Path>,
{
    let mut chroms = Vec::new();
    let mut starts = Vec::new();
    let mut ends = Vec::new();

    for rec in iterator {
        let record = match rec {
            Ok(record) => record,
            Err(VcfParseError::InvariantgVCFLine) => continue,
            Err(err) => return Err(PolarsError::ComputeError(format!("{:?}", err).into())),
        };

        match record.get_span() {
            Ok((start, end)) => {
                chroms.push(record.chrom);
                starts.push(start as i64); // Polars uses i64 for integer columns
                ends.push(end as i64);
            }
            Err(err) => return Err(PolarsError::ComputeError(format!("{:?}", err).into())),
        }
    }

    let mut df = DataFrame::new(vec![
        Series::new("chrom".into(), chroms).into(),
        Series::new("start".into(), starts).into(),
        Series::new("end".into(), ends).into(),
    ])?;

    let file = File::create(output_path)?;
    ParquetWriter::new(file).finish(&mut df)?;
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    let parser = GVcfRecordIterator::from_gzip_path(&args.input)?;
    save_var_regions_as_parquet(parser, &args.output)?;

    println!(
        "Wrote variant spans from {} to {}",
        args.input.display(),
        args.output.display()
    );
    Ok(())
}
