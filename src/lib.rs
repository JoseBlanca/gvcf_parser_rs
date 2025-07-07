use pyo3::prelude::*;
use pyo3::types::PyModule;
use rust_htslib::bgzf::Reader as BgzfReader;
use rust_htslib::tpool::ThreadPool;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Stdin};
use std::path::Path;
use thiserror::Error;

const MISSING_GT: i32 = -1;
const VCF_MIN_COLUMNS: usize = 9;
const CHROM_COLUMN: usize = 0;
const POS_COLUMN: usize = 1;
const REF_ALLELE_COLUMN: usize = 3;
const ALT_ALLELE_COLUMN: usize = 4;
const QUAL_COLUMN: usize = 5;
const FORMAT_COLUMN: usize = 8;
const FIRST_SAMPLE_COLUMN: usize = 9;

#[derive(Error, Debug)]
pub enum VcfParseError {
    #[error("Invalid allele '{allele}'")]
    InvalidAllele { allele: String },

    #[error("Insufficient columns in VCF line: '{line}'")]
    NotEnoughColumns { line: String },

    #[error("Insufficient columns in CHROM header line")]
    NotEnoughColumnsInChromLine,

    #[error("Invalid position value '{value}' in line: '{line}'")]
    InvalidPosition { value: String, line: String },

    #[error("Invalid quality value '{value}': {line}")]
    InvalidQuality { value: String, line: String },

    #[error("Missing GT field in sample '{sample}' in line '{line}'")]
    MissingGtField { sample: String, line: String },

    #[error("FORMAT column (#8) not found in line '{line}'")]
    FormatColumnNotFound { line: String },

    #[error("GT field not found in FORMAT column in line '{line}'")]
    MissingGtFieldInFormat { line: String },

    #[error("Not possible to extract ploidy from line '{line}'")]
    ErrorFindingPloidy { line: String },

    #[error("Inconsistent ploidies found in line '{line}'")]
    InconsistentPloidies { line: String },

    #[error("Observed ({observed}) and given ({given}) ploidies are different line '{line}'")]
    DifferentObservedPloidy {
        line: String,
        observed: usize,
        given: usize,
    },

    #[error("I/O error: {source}")]
    Io {
        #[from]
        source: std::io::Error,
    },

    #[error("I/O error creating the ThreadPool to decompress the VCF file")]
    ThreadPoolError,

    #[error("I/O error opening path: '{path}'")]
    PathError { path: String },

    #[error("Magic byte error")]
    MagicByteError,

    #[error("Gzip in stdin is not supported")]
    GzipInStdinNotSupported,

    #[error("VCF file should be gzipped")]
    VCFFileShouldBeGzipped,
}

pub type VcfResult<T> = std::result::Result<T, VcfParseError>;

fn set_gt(
    genotypes: &mut Vec<i32>,
    sample_idx: usize,
    allele_idx: usize,
    ploidy: usize,
    value: i32,
) {
    let pos = sample_idx * ploidy + allele_idx;
    genotypes[pos] = value;
}

fn parse_allele(allele_str: &str) -> VcfResult<i32> {
    match allele_str {
        "." => Ok(MISSING_GT),
        _ => match allele_str.parse::<i32>() {
            Ok(allele) => Ok(allele),
            Err(_e) => Err(VcfParseError::InvalidAllele {
                allele: (allele_str.to_string()),
            }),
        },
    }
}

fn get_gt_index_from_format_field(cols: &[&str], line: &str) -> VcfResult<usize> {
    cols.get(FORMAT_COLUMN)
        .ok_or(VcfParseError::FormatColumnNotFound {
            line: line.to_string(),
        })?
        .split(":")
        .position(|f| f == "GT")
        .ok_or(VcfParseError::MissingGtFieldInFormat {
            line: line.to_string(),
        })
}

fn look_for_ploidy(line: &str) -> VcfResult<usize> {
    let cols: Vec<&str> = line.trim_end().split('\t').collect();
    let gt_idx = get_gt_index_from_format_field(&cols, line)?;
    for sample_field in &cols[FIRST_SAMPLE_COLUMN..] {
        let gt_str =
            sample_field
                .split(':')
                .nth(gt_idx)
                .ok_or_else(|| VcfParseError::MissingGtField {
                    sample: sample_field.to_string(),
                    line: line.to_string(),
                })?;
        if gt_str == "." {
            continue;
        };
        return Ok(gt_str.split(|c| c == '/' || c == '|').count());
    }
    Err(VcfParseError::ErrorFindingPloidy {
        line: line.to_string(),
    })
}

#[derive(Debug, PartialEq, Eq)]
enum VcfSection {
    Header,
    Body,
}
#[derive(Debug)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub alleles: Vec<String>,
    pub qual: f32,
    pub genotypes: Vec<i32>,
}

impl VcfRecord {
    pub fn from_line(
        num_samples: usize,
        ploidy: usize,
        reference_gt: &str,
        line: &str,
    ) -> VcfResult<Self> {
        let cols: Vec<&str> = line.trim_end().split('\t').collect();
        if cols.len() < VCF_MIN_COLUMNS {
            return Err(VcfParseError::NotEnoughColumns {
                line: (line.to_string()),
            });
        }

        let ref_allele = cols[REF_ALLELE_COLUMN];
        let alt_alleles = cols[ALT_ALLELE_COLUMN];
        let alleles: Vec<String>;
        if alt_alleles == "." {
            alleles = std::iter::once(ref_allele).map(str::to_string).collect();
        } else {
            alleles = std::iter::once(ref_allele)
                .chain(alt_alleles.split(','))
                .map(str::to_string)
                .collect();
        }

        let qual = match cols[QUAL_COLUMN] {
            "." => f32::NAN,
            s => s
                .parse::<f32>()
                .map_err(|_error| VcfParseError::InvalidQuality {
                    value: s.to_string(),
                    line: line.to_string(),
                })?,
        };

        let gt_idx = get_gt_index_from_format_field(&cols, line)?;

        let mut genotypes: Vec<i32> = vec![0; num_samples * ploidy];
        let mut observed_ploidy: Option<usize> = None;
        for (sample_idx, sample_field) in cols[FIRST_SAMPLE_COLUMN..].iter().enumerate() {
            if gt_idx == 0 && sample_field.starts_with(reference_gt) {
                continue;
            };
            let gt_str = sample_field.split(':').nth(gt_idx).ok_or_else(|| {
                VcfParseError::MissingGtField {
                    sample: sample_field.to_string(),
                    line: line.to_string(),
                }
            })?;
            if gt_str == reference_gt {
                observed_ploidy = Some(ploidy);
                continue;
            };

            if gt_str == "." {
                for allele_idx in 0..ploidy {
                    set_gt(&mut genotypes, sample_idx, allele_idx, ploidy, MISSING_GT);
                }
                continue;
            }

            let mut allele_idx: usize = 0;
            for allele_str in gt_str.split(|c| c == '/' || c == '|') {
                if allele_str == "0" {
                    allele_idx += 1;
                    continue;
                };
                set_gt(
                    &mut genotypes,
                    sample_idx,
                    allele_idx,
                    ploidy,
                    parse_allele(allele_str)?,
                );
                allele_idx += 1;
            }
            if let Some(value) = observed_ploidy {
                if value != allele_idx {
                    return Err(VcfParseError::InconsistentPloidies {
                        line: line.to_string(),
                    });
                }
            } else {
                observed_ploidy = Some(allele_idx);
            }
        }

        if let Some(value) = observed_ploidy {
            if value != ploidy {
                return Err(VcfParseError::DifferentObservedPloidy {
                    line: line.to_string(),
                    observed: value,
                    given: ploidy,
                });
            }
        }

        Ok(VcfRecord {
            chrom: cols[CHROM_COLUMN].to_string(),
            pos: cols[POS_COLUMN]
                .parse()
                .map_err(|_e| VcfParseError::InvalidPosition {
                    value: cols[POS_COLUMN].to_string(),
                    line: line.to_string(),
                })?,
            alleles,
            qual,
            genotypes,
        })
    }
}

pub struct VcfRecordIterator<R: BufRead> {
    reader: R,
    line: String,
    section: VcfSection,
    num_samples: usize,
    ploidy: usize,
    reference_gt: String,
}

impl<R: BufRead> VcfRecordIterator<R> {
    fn new(reader: R) -> Self {
        VcfRecordIterator {
            reader,
            line: String::new(),
            section: VcfSection::Header,
            num_samples: 0,
            ploidy: 0,
            reference_gt: String::new(),
        }
    }

    fn parse_variant(&self) -> VcfResult<VcfRecord> {
        VcfRecord::from_line(
            self.num_samples,
            self.ploidy,
            &self.reference_gt,
            &self.line,
        )
    }

    fn process_chrom_line(&mut self) -> Option<VcfResult<VcfRecord>> {
        let fields_: Vec<&str> = self.line.trim_end().split('\t').collect();
        if fields_.len() < VCF_MIN_COLUMNS {
            return Some(Err(VcfParseError::NotEnoughColumnsInChromLine));
        }
        self.num_samples = fields_.len() - VCF_MIN_COLUMNS;
        None // Continue processing, don't return a record yet
    }

    fn process_first_variant_line(&mut self) -> Option<VcfResult<VcfRecord>> {
        // Set up ploidy and reference genotype
        match look_for_ploidy(&self.line) {
            Ok(ploidy) => self.ploidy = ploidy,
            Err(e) => return Some(Err(e)),
        }

        self.reference_gt = vec!["0"; self.ploidy].join("/");
        self.section = VcfSection::Body;

        // Parse and return the first variant record
        let record = VcfRecord::from_line(
            self.num_samples,
            self.ploidy,
            &self.reference_gt,
            &self.line,
        );
        self.section = VcfSection::Body;
        Some(record)
    }

    fn process_header_and_first_variant(&mut self) -> Option<VcfResult<VcfRecord>> {
        loop {
            match () {
                _ if self.line.starts_with("##") => None, // Continue
                _ if self.line.starts_with("#CHROM") => self.process_chrom_line(),
                _ => {
                    return self.process_first_variant_line();
                }
            };

            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    continue;
                }
                Err(e) => return Some(Err(VcfParseError::Io { source: e })),
            }
        }
    }
}

impl<R: BufRead> Iterator for VcfRecordIterator<R> {
    type Item = VcfResult<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line.clear();

        let result = match self.reader.read_line(&mut self.line) {
            Ok(0) => None, // EOF
            Ok(_) => match self.section {
                VcfSection::Body => return Some(self.parse_variant()),
                VcfSection::Header => return self.process_header_and_first_variant(),
            },
            Err(error) => Some(Err(VcfParseError::from(error))),
        };
        return result;
    }
}

impl<R: Read> VcfRecordIterator<BufReader<R>> {
    pub fn from_reader(reader: R) -> Self {
        let buf_reader = BufReader::new(reader);
        VcfRecordIterator::new(buf_reader)
    }
}

impl VcfRecordIterator<BufReader<rust_htslib::bgzf::Reader>> {
    pub fn from_gzipped_vcf_path<P: AsRef<Path>>(
        path: P,
        n_threads: u32,
    ) -> VcfResult<(Self, Option<ThreadPool>)> {
        let file = File::open(&path)?;
        let mut buf_reader = BufReader::new(file);

        let num_bytes = 4;
        let buffer = buf_reader.fill_buf()?;
        let first_bytes = &buffer[..num_bytes.min(buffer.len())];
        let file_is_gzziped =
            are_gzipped_magic_bytes(first_bytes).map_err(|_| VcfParseError::MagicByteError)?;
        if !file_is_gzziped {
            return Err(VcfParseError::VCFFileShouldBeGzipped);
        }

        let mut bgz_reader =
            BgzfReader::from_path(&path).map_err(|_e| VcfParseError::PathError {
                path: path.as_ref().to_string_lossy().into_owned(),
            })?;
        let pool = ThreadPool::new(n_threads).map_err(|_e| VcfParseError::ThreadPoolError)?;
        bgz_reader
            .set_thread_pool(&pool)
            .map_err(|_e| VcfParseError::ThreadPoolError)?;
        let buf_bgz_reader = BufReader::new(bgz_reader);
        let parser = VcfRecordIterator::new(buf_bgz_reader);
        Ok((parser, Some(pool)))
    }
}

impl VcfRecordIterator<BufReader<std::io::Stdin>> {
    pub fn from_stdin() -> VcfResult<VcfRecordIterator<BufReader<Stdin>>> {
        let stdin = std::io::stdin();
        let mut reader = BufReader::new(stdin);

        let buffer = reader.fill_buf()?;
        let first_bytes = &buffer[..4.min(buffer.len())];
        let file_is_gzziped =
            are_gzipped_magic_bytes(first_bytes).map_err(|_| VcfParseError::MagicByteError)?;
        if file_is_gzziped {
            return Err(VcfParseError::GzipInStdinNotSupported);
        }
        Ok(Self::new(reader))
    }
}

#[derive(Error, Debug)]
pub enum MagicByteError {
    #[error("Insufficient bytes: got {got}, need at least {need}")]
    InsufficientBytes { got: usize, need: usize },
}

fn are_gzipped_magic_bytes(first_bytes: &[u8]) -> Result<bool, MagicByteError> {
    if first_bytes.len() < 2 {
        return Err(MagicByteError::InsufficientBytes {
            got: first_bytes.len(),
            need: 2,
        });
    }
    Ok(first_bytes[0] == 0x1f && first_bytes[1] == 0x8b)
}

#[pyclass(name = "VcfRecord")]
#[derive(Debug, Clone)]
pub struct PyVcfRecord {
    #[pyo3(get)]
    pub chrom: String,
    #[pyo3(get)]
    pub pos: u32,
    #[pyo3(get)]
    pub alleles: Vec<String>,
    #[pyo3(get)]
    pub qual: f32,
    #[pyo3(get)]
    pub genotypes: Vec<i32>,
}

impl From<VcfRecord> for PyVcfRecord {
    fn from(rec: VcfRecord) -> Self {
        PyVcfRecord {
            chrom: rec.chrom,
            pos: rec.pos,
            alleles: rec.alleles,
            qual: rec.qual,
            genotypes: rec.genotypes,
        }
    }
}

#[pyclass(name = "VcfRecordIterator", unsendable)]
pub struct PyVcfRecordIterator {
    inner: Box<dyn Iterator<Item = VcfResult<VcfRecord>>>,
    _pool: Option<ThreadPool>,
}

#[pymethods]
impl PyVcfRecordIterator {
    #[new]
    fn new(path: String, n_threads: u32) -> PyResult<Self> {
        let (parser, pool) = VcfRecordIterator::from_gzipped_vcf_path(path, n_threads)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("{:?}", e)))?;
        Ok(Self {
            inner: Box::new(parser),
            _pool: pool,
        })
    }

    fn __iter__(slf: PyRefMut<'_, Self>) -> Py<PyVcfRecordIterator> {
        slf.into()
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyResult<PyVcfRecord>> {
        slf.inner.next().map(|result| match result {
            Ok(rec) => Ok(PyVcfRecord::from(rec)),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "{:?}",
                e
            ))),
        })
    }
}

#[pymodule]
fn vcfparser(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyVcfRecordIterator>()?;
    m.add_class::<PyVcfRecord>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    const SAMPLE_VCF: &str = "##fileformat=VCFv4.5
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">
##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/6000:43:5:.,.
20\t17330\t.\tT\tA\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t.:35:4\t0/2:17:2\t./1:40:3";

    #[test]
    //#[ignore]
    fn test_parse_vcf_iter() {
        let reader = BufReader::new(SAMPLE_VCF.as_bytes());
        let parser = VcfRecordIterator::from_reader(reader);

        let mut count = 0;
        for record_result in parser {
            match record_result {
                Ok(_record) => {
                    count += 1;
                }
                Err(e) => {
                    eprintln!("Error parsing record {}: {}", count, e);
                    break;
                }
            }
        }
        assert_eq!(count, 6);
    }

    #[test]
    //#[ignore]
    fn test_parse_vcf_gz_file_iter() -> Result<(), Box<dyn std::error::Error>> {
        let result = VcfRecordIterator::from_gzipped_vcf_path(
            "/home/jose/analyses/g2psol/source_data/TS.vcf.gz",
            4,
        )?;
        let (parser, _pool) = result;

        let mut count = 0;
        for record_result in parser {
            match record_result {
                Ok(_record) => {
                    count += 1;
                }
                Err(e) => {
                    eprintln!("Error parsing record {}: {}", count, e);
                    break;
                }
            }
        }
        assert_eq!(count, 4249);
        Ok(())
    }
}
