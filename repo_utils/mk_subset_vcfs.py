
import pysam

infn="/Users/Mac/Downloads/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
v = pysam.VariantFile(infn)

for chunk in range(3):
    o = pysam.VariantFile(f"chunk{chunk}.vcf", 'w', header=v.header)
    for i in range(1000):
        o.write(next(v))
