import sys
import pysam
import numpy as np
import pysam
import truvari
import logging
import allel
import joblib


truvari.setup_logging()
logging.info("Starting")
fn = sys.argv[1]
# fields = optionaly AF (make sure its consistent with Number=[./1 whatever]
d = joblib.dump(allel.read_vcf(fn, fields=["calldata/GT", "samples"]) , "test.jl")

logging.info("Read")
import h5py

#d = h5py.File("test.hdf5")
d = joblib.load("test.jl")
vcf_samples = d["samples"]
gts = allel.GenotypeArray(d["calldata/GT"])

is_het = gts.is_het()
is_hom = gts.is_hom_alt()
v_count = is_het | is_hom
# topN approach - choose the sample with the most
# Then choose the sample with the most without this sample

max_reporting = 100
num_samples = gts.shape[1]
num_vars = gts.shape[0]

variant_mask = np.zeros(num_vars, dtype='bool') # used variants
sample_mask = np.zeros(num_samples, dtype='bool') # used samples

logging.info(f"Sample Count {num_samples}")
logging.info(f"Variant Count {num_vars}")

# Okay, I accidentally implemented the greedy approach... 
# So I'm curious how topN works, but until then, whatever
print("sample\tvarcount\tnew_count\ttot_captured\tpct_captured")
for i in range(max_reporting): # picking the best N samples
    # of the variants remaining
    cur_view = v_count[~variant_mask]
    # how many variants per sample
    cur_sample_count = cur_view.sum(axis=0)
    # use the sample with the most variants (excluding used samples)
    use_sample = np.argmax(cur_sample_count)

    # how many does this sample have overall?
    use_sample_variant_count = v_count[:,use_sample].sum()

    # number of new variants added
    new_variant_count = np.max(cur_sample_count[~sample_mask])
    # don't want to use these variants anymore
    variant_mask = variant_mask | v_count[:, use_sample]
    # or this sample
    sample_mask[use_sample] = True

    # our running total number of variants
    upto_now = variant_mask.sum()
    print(vcf_samples[use_sample], 
          use_sample_variant_count,
          new_variant_count,
          upto_now,
          round(upto_now / num_vars, 4), sep='\t')

logging.info("Finished")
