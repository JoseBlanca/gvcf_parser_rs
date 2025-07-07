import vcfparser
from tempfile import NamedTemporaryFile, TemporaryDirectory
from subprocess import run

VCF1 = """##fileformat=VCFv4.5
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t10\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/6000:43:5:.,.
20\t20\t.\tT\tA\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t30\t.\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t40\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t50\tindel\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
"""
VCF2 = """##fileformat=VCFv4.5
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00004\tNA00005\tNA00006
20\t10\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/6000:43:5:.,.
20\t20\t.\tT\tA\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t30\t.\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t40\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t51\tindel\tT\tA,C\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
"""


def create_temp_vcf(content, prefix, temp_dir):
    vcf = NamedTemporaryFile(suffix=".vcf", prefix=prefix, dir=temp_dir)
    with open(vcf.name, "wt") as fhand:
        fhand.write(content)
        fhand.flush()
    cmd = ["bgzip", vcf.name]
    run(cmd, check=True, capture_output=True)
    return vcf


def setup_vcfs():
    temp_dir = TemporaryDirectory(dir=".")
    temp_vcf1 = create_temp_vcf(VCF1, "vcf1", temp_dir=temp_dir.name)
    temp_vcf2 = create_temp_vcf(VCF2, "vcf2", temp_dir=temp_dir.name)
    return temp_vcf1, temp_vcf2, temp_dir


def calc_var_len(var):
    return max(len(allele) for allele in var.alleles)


def calc_var_span(var):
    var_len = calc_var_len(var)
    return var.pos, var.pos + var_len - 1


def main():
    temp_vcf1, temp_vcf2, temp_dir = setup_vcfs()
    vars1 = vcfparser.PyVcfRecordIterator(temp_vcf1.name + ".gz", 1)
    vars2 = vcfparser.PyVcfRecordIterator(temp_vcf2.name + ".gz", 1)
    for var in vars1:
        print(calc_var_span(var))


main()
