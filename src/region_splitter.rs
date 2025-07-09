use crate::errors::VcfParseError;
use crate::utils_magic::file_is_gzipped;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub type VcfResult<T> = std::result::Result<T, VcfParseError>;

const NON_REF: &str = "<NON_REF>";

pub struct GVcfRecord {
    pub chrom: String,
    pub pos: u32,
    alleles: Vec<String>,
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
}

impl<B: BufRead> GVcfRecordIterator<B> {
    fn new(reader: B) -> Self {
        GVcfRecordIterator {
            reader: reader,
            line: String::new(),
            section: VcfSection::Header,
        }
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

impl<R: BufRead> Iterator for GVcfRecordIterator<R> {
    type Item = VcfResult<GVcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line.clear();

        match self.reader.read_line(&mut self.line) {
            Ok(0) => return None, // EOF
            Ok(_) => {
                if self.section == VcfSection::Header {
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
                    if !self.line.starts_with("#CHROM") {
                        return Some(Err(VcfParseError::MalformedHeader));
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
                } else {
                    Some(GVcfRecord::from_line(&self.line))
                }
            }
            Err(error) => Some(Err(VcfParseError::from(error))),
        }
    }
}
