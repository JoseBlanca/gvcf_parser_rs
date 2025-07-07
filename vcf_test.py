import vcfparser
from cyvcf2 import VCF
import gzip

fpath = "/home/jose/analyses/g2psol/source_data/TS.vcf.gz"
if True:
    vcf_iterator = vcfparser.PyVcfRecordIterator(fpath, 4)
    print(sum([1 for record in vcf_iterator]))
    for record in vcfparser.PyVcfRecordIterator(fpath, 4):
        print(record.chrom, record.pos, record.alleles)

elif False:
    print(sum([1 for variant in VCF(fpath, threads=4)]))
else:
    print(sum(1 for line in gzip.open(fpath, "rt")))
