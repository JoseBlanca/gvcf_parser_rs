use std::error::Error;
use std::io::BufRead;

const MISSING_GT: i32 = -1;

#[derive(Debug)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub alleles: Vec<String>,
    pub qual: f32,
    pub genotypes: Vec<i32>,
}

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

fn digit_str_to_int(s: &str) -> i32 {
    match s {
        "." => MISSING_GT,
        "0" => 0,
        "1" => 1,
        "2" => 2,
        "3" => 3,
        "4" => 4,
        "5" => 5,
        "6" => 6,
        "7" => 7,
        "8" => 8,
        "9" => 9,
        allele => allele.parse::<i32>().unwrap_or(MISSING_GT),
    }
}

impl VcfRecord {
    pub fn from_line(
        num_samples: &usize,
        ploidy: usize,
        reference_gt: &str,
        line: &str,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let cols: Vec<&str> = line.trim_end().split('\t').collect();
        if cols.len() < 8 {
            return Err("Not enough columns in VCF line".into());
        }

        let ref_allele = cols[3];
        let alt_alleles = cols[4];
        let alleles: Vec<String>;
        if alt_alleles == "." {
            alleles = std::iter::once(ref_allele).map(str::to_string).collect();
        } else {
            alleles = std::iter::once(ref_allele)
                .chain(alt_alleles.split(','))
                .map(str::to_string)
                .collect();
        }

        let qual = match cols[5] {
            "." => f32::NAN,
            s => s.parse::<f32>()?,
        };

        let gt_idx = get_gt_index_from_format_field(&cols)?;

        let mut genotypes: Vec<i32> = vec![0; num_samples * ploidy];
        let mut observed_ploidy: Option<usize> = None;
        for (sample_idx, sample_field) in cols[9..].iter().enumerate() {
            if gt_idx == 0 && sample_field.starts_with(reference_gt) {
                continue;
            };
            let gt_str = sample_field.split(':').nth(gt_idx).ok_or_else(|| {
                format!(
                    "Missing GT field at index {} in sample '{}'",
                    gt_idx, sample_field
                )
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
                    digit_str_to_int(allele_str),
                );
                allele_idx += 1;
            }
            if let Some(value) = observed_ploidy {
                if value != allele_idx {
                    return Err("Inconsistent observed ploidies in a variant".into());
                }
            } else {
                observed_ploidy = Some(allele_idx);
            }
        }

        if let Some(value) = observed_ploidy {
            if value != ploidy {
                return Err("observed ploidy different than given ploidy".into());
            }
        }

        Ok(VcfRecord {
            chrom: cols[0].to_string(),
            pos: cols[1].parse()?,
            alleles,
            qual,
            genotypes,
        })
    }
}

fn get_gt_index_from_format_field(cols: &[&str]) -> Result<usize, Box<dyn Error>> {
    cols.get(8)
        .ok_or("FORMAT column (#8) not found")?
        .split(":")
        .position(|f| f == "GT")
        .ok_or("GT field not found in FORMAT".into())
}

fn look_for_ploidy(line: &str) -> Result<usize, Box<dyn Error>> {
    let cols: Vec<&str> = line.trim_end().split('\t').collect();
    let gt_idx = get_gt_index_from_format_field(&cols)?;
    for sample_field in &cols[9..] {
        let gt_str = sample_field.split(':').nth(gt_idx).ok_or_else(|| {
            format!(
                "Missing GT field at index {} in sample '{}'",
                gt_idx, sample_field
            )
        })?;
        if gt_str == "." {
            continue;
        };
        return Ok(gt_str.split(|c| c == '/' || c == '|').count());
    }
    Err("It was not possible to determine ploidy".into())
}

#[derive(Debug, PartialEq, Eq)]
enum VcfSection {
    Header,
    Body,
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
    pub fn new(reader: R) -> Self {
        VcfRecordIterator {
            reader,
            line: String::new(),
            section: VcfSection::Header,
            num_samples: 0,
            ploidy: 0,
            reference_gt: String::new(),
        }
    }
}

impl<R: BufRead> Iterator for VcfRecordIterator<R> {
    type Item = Result<VcfRecord, Box<dyn Error>>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line.clear();

        match self.reader.read_line(&mut self.line) {
            Ok(0) => None, // EOF
            Ok(_) => {
                match self.section {
                    VcfSection::Body => {
                        match VcfRecord::from_line(
                            &self.num_samples,
                            self.ploidy,
                            &self.reference_gt,
                            &self.line,
                        ) {
                            Ok(record) => return Some(Ok(record)),
                            Err(e) => return Some(Err(e)),
                        }
                    }
                    VcfSection::Header => {
                        loop {
                            if self.line.starts_with("##") {
                                // Skip header lines, continue to next iteration
                            } else if self.line.starts_with("#CHROM") {
                                let fields_: Vec<&str> = self.line.trim_end().split('\t').collect();
                                if fields_.len() < 9 {
                                    return Some(Err(
                                        "Not enough fields found in the CHROM line".into()
                                    ));
                                }
                                self.num_samples = fields_.len() - 9;
                            } else {
                                match look_for_ploidy(&self.line) {
                                    Ok(ploidy) => self.ploidy = ploidy,
                                    Err(e) => return Some(Err(e)),
                                };
                                self.reference_gt = vec!["0"; self.ploidy].join("/");
                                self.section = VcfSection::Body;
                                let record = VcfRecord::from_line(
                                    &self.num_samples,
                                    self.ploidy,
                                    &self.reference_gt,
                                    &self.line,
                                );
                                self.line.clear();
                                return Some(record);
                            }
                            self.line.clear();
                            match self.reader.read_line(&mut self.line) {
                                Ok(0) => return None, // EOF
                                Ok(_) => {
                                    continue;
                                }
                                Err(e) => return Some(Err(Box::new(e))),
                            }
                        }
                    }
                };
            }
            Err(e) => Some(Err(Box::new(e))),
        }
    }
}

pub fn parse_vcf_iter<R: BufRead>(reader: R) -> VcfRecordIterator<R> {
    VcfRecordIterator::new(reader)
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
        let parser = parse_vcf_iter(reader);

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
        use rust_htslib::bgzf::Reader as BgzfReader;
        use rust_htslib::tpool::ThreadPool;
        let n_threads = 4;
        let pool = ThreadPool::new(n_threads)?;
        let file_name = "/home/jose/analyses/g2psol/source_data/TS.vcf.gz";
        let mut raw_reader = BgzfReader::from_path(file_name)?;
        raw_reader.set_thread_pool(&pool)?;

        let reader = BufReader::new(raw_reader);
        let parser = parse_vcf_iter(reader);

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
