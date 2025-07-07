use rust_htslib::bgzf::Reader as BgzfReader;
use rust_htslib::tpool::ThreadPool;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Stdin};
use std::path::Path;

use crate::errors::VcfParseError;
use crate::utils_magic::are_gzipped_magic_bytes;

const MISSING_GT: i32 = -1;
const VCF_MIN_COLUMNS: usize = 9;
const CHROM_COLUMN: usize = 0;
const POS_COLUMN: usize = 1;
const REF_ALLELE_COLUMN: usize = 3;
const ALT_ALLELE_COLUMN: usize = 4;
const QUAL_COLUMN: usize = 5;
const FORMAT_COLUMN: usize = 8;
const FIRST_SAMPLE_COLUMN: usize = 9;

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
