use std::error::Error;
use std::io::BufRead;

#[derive(Debug)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    pub qual: f32,
    pub samples: Vec<String>,
}

pub fn parse_vcf<R: BufRead>(reader: R) -> Result<Vec<VcfRecord>, Box<dyn Error>> {
    for (line_idx, line_result) in reader.lines().enumerate() {
        match line_result {
            Ok(line) => println!(
                "OK {line_idx} len={}: {}",
                line.len(),
                &line.chars().take(20).collect::<String>()
            ),
            Err(e) => {
                println!("Error at line {line_idx}: {e}");
                return Err(Box::new(e));
            }
        }
        //println!("{line_idx}");
    }
    Err(format!("Hola").into())
}

pub fn parse_vcf2<R: BufRead>(mut reader: R) -> Result<Vec<VcfRecord>, Box<dyn Error>> {
    let mut line = String::new();
    let mut records: Vec<VcfRecord> = Vec::new();

    let mut idx = 0;
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }

        if line.starts_with("##") {
            continue;
        }
        println!(
            "OK {idx} len={}: {}",
            line.len(),
            &line.chars().take(20).collect::<String>()
        );
        idx += 1;
    }

    Err("Finished reading".into())
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
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
20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3";

    #[test]
    #[ignore]
    fn test_parse_vcf() {
        let reader = BufReader::new(SAMPLE_VCF.as_bytes());
        let _res = parse_vcf(reader);
    }
    #[test]
    fn test_parse_vcf_gz_file() {
        let file_name = "/home/jose/analyses/g2psol/source_data/TS.vcf.gz";
        let file_name = "/home/jose/devel/vcf_parser_rs/all.vcf.gz";

        let file = File::open(file_name).expect("Failed to open VCF file");
        let gz = flate2::read::GzDecoder::new(file);
        let reader = BufReader::new(gz);
        let _res = parse_vcf2(reader);
    }
}
