use crate::errors::VcfParseError;
use crate::utils_magic::{file_is_bgzipped, file_is_gzipped};
use flate2::read::MultiGzDecoder;
use rust_htslib::bgzf::Reader as BgzfReader;
use rust_htslib::tpool::ThreadPool;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub type VcfResult<T> = std::result::Result<T, VcfParseError>;

const NON_REF: &str = "<NON_REF>";

#[derive(Debug)]
pub struct GVcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub alleles: Vec<String>,
}

impl GVcfRecord {
    fn from_line(line: &str) -> VcfResult<Self> {
        let mut fields = line.splitn(6, '\t');
        let chrom = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;

        let pos = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;
        fields.next();
        let ref_allele = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;
        let alt_alleles = fields
            .next()
            .ok_or_else(|| VcfParseError::GVCFLineNotEnoughFields)?;

        if alt_alleles == NON_REF {
            return Err(VcfParseError::InvariantgVCFLine);
        }

        let pos = pos
            .parse::<u32>()
            .map_err(|_| VcfParseError::GVCFLineNotEnoughFields)?;

        let alleles: Vec<String> = std::iter::once(ref_allele)
            .chain(alt_alleles.split(','))
            .filter(|allele| allele != &NON_REF)
            .map(str::to_string)
            .collect();

        Ok(GVcfRecord {
            chrom: chrom.to_string(),
            pos: pos,
            alleles: alleles,
        })
    }
    pub fn get_span(self: &GVcfRecord) -> VcfResult<(u32, u32)> {
        let max_allele_len = self.alleles.iter().map(|allele| allele.len()).max().ok_or(
            VcfParseError::RuntimeError {
                message: "There should be at least one allele".to_string(),
            },
        )?;
        if max_allele_len == 1 {
            Ok((self.pos, self.pos))
        } else {
            Ok((self.pos, self.pos + max_allele_len as u32 - 1))
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
enum VcfSection {
    Header,
    Body,
}

pub struct GVcfRecordIterator<B: BufRead> {
    reader: B,
    line: String,
    section: VcfSection,
    buffer: VecDeque<GVcfRecord>,
}

impl<B: BufRead> GVcfRecordIterator<B> {
    fn new(reader: B) -> Self {
        GVcfRecordIterator {
            reader: reader,
            line: String::new(),
            section: VcfSection::Header,
            buffer: VecDeque::new(),
        }
    }
    fn process_header_and_first_variant(&mut self) -> Option<VcfResult<GVcfRecord>> {
        loop {
            if self.line.starts_with("##") {
                self.line.clear();
                match self.reader.read_line(&mut self.line) {
                    Ok(0) => return Some(Err(VcfParseError::BrokenHeader)),
                    Ok(_) => {
                        if !self.line.starts_with("##") {
                            break;
                        }
                    }
                    Err(error) => return Some(Err(VcfParseError::from(error))),
                }
            }
        }
        self.line.clear();
        match self.reader.read_line(&mut self.line) {
            Ok(0) => None, // EOF
            Ok(_) => {
                self.section = VcfSection::Body;
                Some(GVcfRecord::from_line(&self.line))
            }
            Err(error) => return Some(Err(VcfParseError::from(error))),
        }
    }
    pub fn fill_buffer(&mut self, n_items: usize) -> VcfResult<usize> {
        let mut n_items_added: usize = 0;
        while self.buffer.len() < n_items {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    if self.section == VcfSection::Header {
                        let result = self.process_header_and_first_variant();
                        if let Some(Ok(record)) = result {
                            self.buffer.push_back(record);
                            n_items_added += 1;
                        }
                    } else {
                        match GVcfRecord::from_line(&self.line) {
                            Ok(record) => {
                                self.buffer.push_back(record);
                                n_items_added += 1;
                            }
                            Err(VcfParseError::InvariantgVCFLine) => continue, // skip
                            Err(err) => return Err(err),
                        }
                    }
                }
                Err(err) => {
                    return Err(VcfParseError::from(err));
                }
            }
        }
        Ok(n_items_added)
    }

    pub fn peek_items_in_buffer(&self) -> impl Iterator<Item = &GVcfRecord> {
        self.buffer.iter()
    }
}

impl<R: Read> GVcfRecordIterator<BufReader<R>> {
    pub fn from_reader(reader: R) -> Self {
        let buf_reader = BufReader::new(reader);
        GVcfRecordIterator::new(buf_reader)
    }
}
impl<R: Read> GVcfRecordIterator<BufReader<MultiGzDecoder<R>>> {
    pub fn from_gzip_reader(reader: R) -> Self {
        let gz_decoder = MultiGzDecoder::new(reader);
        let buf_reader = BufReader::new(gz_decoder);
        GVcfRecordIterator::new(buf_reader)
    }
}
impl GVcfRecordIterator<BufReader<MultiGzDecoder<File>>> {
    pub fn from_gzip_path<P: AsRef<Path>>(path: P) -> VcfResult<Self> {
        if !file_is_gzipped(&path).map_err(|_| VcfParseError::MagicByteError)? {
            return Err(VcfParseError::VCFFileShouldBeGzipped);
        }
        let file = File::open(&path)?;
        let gz_decoder = MultiGzDecoder::new(file);
        let buf_reader = BufReader::new(gz_decoder);
        Ok(GVcfRecordIterator::new(buf_reader))
    }
}

impl GVcfRecordIterator<BufReader<rust_htslib::bgzf::Reader>> {
    pub fn from_bgzip_path<P: AsRef<Path>>(
        path: P,
        n_threads: u32,
    ) -> VcfResult<(Self, ThreadPool)> {
        if !file_is_bgzipped(&path).map_err(|_| VcfParseError::MagicByteError)? {
            return Err(VcfParseError::VCFFileShouldBeBGzipped);
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
        Ok((GVcfRecordIterator::new(buf_bgz_reader), pool))
    }
}

impl<R: BufRead> Iterator for GVcfRecordIterator<R> {
    type Item = VcfResult<GVcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line.clear();

        match self.reader.read_line(&mut self.line) {
            Ok(0) => return None, // EOF
            Ok(_) => match self.section {
                VcfSection::Body => Some(GVcfRecord::from_line(&self.line)),
                VcfSection::Header => self.process_header_and_first_variant(),
            },
            Err(error) => Some(Err(VcfParseError::from(error))),
        }
    }
}

// fn get_span_covers_at_least(chrom, end) -> n_records, goes_beyond:bool, new_end
