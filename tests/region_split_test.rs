use gvcfparser::{
    errors::VcfParseError,
    gvcf_parser::{GVcfRecord, GVcfRecordIterator},
};
use std::fs::File;
use std::io::BufReader;

const SAMPLE_GVCF: &str = "##
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t14370\t.\tG\t<NON_REF>\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/6000:43:5:.,.
20\t17330\t.\tT\tA,<NON_REF>\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t17331\t.\tA\tG,T,<NON_REF>\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t17332\t.\tT\t<NON_REF>\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t17333\t.\tGTC\tG,GTCT,<NON_REF>\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
20\t17334\t.\tGTC\tG,GTCT,<NON_REF>\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t.:35:4\t0/2:17:2\t./1:40:3";

#[test]
fn test_gvcf_parsing() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let parser = GVcfRecordIterator::from_reader(reader);
    let mut n_variants: u32 = 0;
    for record in parser {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    assert_eq!(n_variants, 4);
}

#[test]
fn test_buffer() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let mut var_iterator = GVcfRecordIterator::from_reader(reader);
    assert!(matches!(var_iterator.fill_buffer(3), Ok(3)));
    assert!(matches!(var_iterator.fill_buffer(1), Ok(0)));
    assert!(matches!(var_iterator.fill_buffer(4), Ok(1)));
    assert!(matches!(var_iterator.fill_buffer(5), Ok(0)));
    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17330);
    let buffered_items = var_iterator.peek_items_in_buffer();
    let poss = [17331, 17333, 17334];
    for (expected_pos, variant) in poss.iter().zip(buffered_items) {
        assert_eq!(&variant.pos, expected_pos);
    }
}
#[test]
fn test_buffer2() {
    let reader = BufReader::new(SAMPLE_GVCF.as_bytes());
    let mut var_iterator = GVcfRecordIterator::from_reader(reader);
    assert!(matches!(var_iterator.fill_buffer(2), Ok(2)));
    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17330);
    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17331);

    let buffered_items: Vec<&GVcfRecord> = var_iterator.peek_items_in_buffer().collect();
    assert_eq!(buffered_items.len(), 0);

    let variant = var_iterator.next().unwrap().unwrap();
    assert_eq!(variant.pos, 17333);

    assert!(matches!(var_iterator.fill_buffer(2), Ok(0)));

    let buffered_items = var_iterator.peek_items_in_buffer();
    let poss = [17334];
    for (expected_pos, variant) in poss.iter().zip(buffered_items) {
        assert_eq!(&variant.pos, expected_pos);
    }
}

#[test]
fn test_gzip_reader() {
    let file = File::open("tests/data/sample.g.vcf.gz").expect("Problem opening test file");
    let records = GVcfRecordIterator::from_gzip_reader(file);

    let mut n_variants: u32 = 0;
    for record in records {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    assert_eq!(n_variants, 0);
}
#[test]
fn test_gzip_path() {
    let path = "tests/data/sample.g.vcf.gz";
    let records = GVcfRecordIterator::from_gzip_path(path).expect("Problem opening test file");

    let mut n_variants: u32 = 0;
    for record in records {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    assert_eq!(n_variants, 0);
}

#[test]
fn test_performance() {
    let path = "sample_files/sample.g.vcf.gz";
    //let records = GVcfRecordIterator::from_gzip_path(path).expect("Problem opening test file");
    //println!("g.vcf.gzip");

    let path = "sample_files/sample.g.vcf.bgz";
    let records = GVcfRecordIterator::from_gzip_path(path).expect("Problem opening test file");
    println!("{}", path);

    let mut n_variants: u32 = 0;
    let mut n_invariants: u32 = 0;
    for record in records {
        match record {
            Ok(_variant) => {
                n_variants += 1;
            }
            Err(VcfParseError::InvariantgVCFLine) => {
                n_invariants += 1;
            }
            Err(error) => {
                //Fail test
                panic!("Unexpected error: {}", error);
            }
        }
    }
    println!("Num. variant loci: {n_variants}");
    println!("Num. invariant loci: {n_invariants}");
    println!("Num.loci: {}", n_invariants + n_variants);
}

#[test]
fn test_g_vcf_record() {
    let pos: u32 = 10;
    let alleles = vec!["A".to_string(), "C".to_string()];
    let snp = GVcfRecord {
        chrom: "chr1".to_string(),
        pos: pos,
        alleles: alleles,
    };
    assert!(matches!(snp.get_span(), Ok((10, 10))));

    let alleles = vec!["AT".to_string(), "A".to_string()];
    let snp = GVcfRecord {
        chrom: "chr1".to_string(),
        pos: pos,
        alleles: alleles,
    };
    assert!(matches!(snp.get_span(), Ok((10, 11))));

    let alleles = vec!["A".to_string(), "ATT".to_string()];
    let snp = GVcfRecord {
        chrom: "chr1".to_string(),
        pos: pos,
        alleles: alleles,
    };
    assert!(matches!(snp.get_span(), Ok((10, 12))));
}
