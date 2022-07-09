
import pysam

v = pysam.VariantFile("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

for chunk in range(2):
    o = pysam.VariantFile(f"chunk{chunk}.vcf", 'w', header=v.header)
    for i in range(1000):
        o.write(next(v))
