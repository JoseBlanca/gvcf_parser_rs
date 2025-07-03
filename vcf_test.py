import vcfparser
from cyvcf2 import VCF
from pigz import PigzFile

fpath = "/home/jose/analyses/g2psol/source_data/TS.vcf.gz"
if False:
    vcf_iterator = vcfparser.PyVcfRecordIterator(fpath, 4)
    print(sum([1 for record in vcf_iterator]))
elif False:
    print(sum([1 for variant in VCF(fpath, threads=4)]))
