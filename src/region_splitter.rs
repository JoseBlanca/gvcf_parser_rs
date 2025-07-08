use std::io::{BufRead, BufReader, Read};

use crate::errors::VcfParseError;

pub type VcfResult<T> = std::result::Result<T, VcfParseError>;

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

        if alt_alleles == "<NO_REF>" {
            return Err(VcfParseError::InvariantgVCFLine);
        }

        let pos = pos
            .parse::<u32>()
            .map_err(|_| VcfParseError::GVCFLineNotEnoughFields)?;

        let alleles: Vec<String> = std::iter::once(ref_allele)
            .chain(alt_alleles.split(','))
            .filter(|allele| allele != &"<NO_REF>")
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
