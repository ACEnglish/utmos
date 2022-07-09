import logging
import argparse

import joblib
import truvari
import numpy as np

from utmos.convert import read_vcf

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="select", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_files", nargs="+", type=str,
                        help="Input VCF or jl files")
    parser.add_argument("-o", "--out", type=str, default="/dev/stdout",
                        help="Output file (stdout)")
    parser.add_argument("-k", "--keep", type=float, default=0.02,
                        help="Number of samples to select as a percent if <1 or count if >=1 (%(default)s)")
    parser.add_argument("--safe", action='store_true',
                        help="Ensure input files have same sample names")
    args = parser.parse_args(args)
    truvari.setup_logging()
    return args


def calculate(data, out_fn, max_reporting=0.02):
    """
    Do the selection calculation
    if max_reporting [0,1], select that percent of samples
    if max_reporting >= 1, select that number of samples
    """
    v_count = data['GT']
    vcf_samples = data['samples']

    num_samples = v_count.shape[1]
    num_vars = v_count.shape[0]

    if max_reporting < 1:
        max_reporting = int(num_samples * max_reporting)
    else:
        max_reporting = int(max_reporting)

    variant_mask = np.zeros(num_vars, dtype='bool') # used variants
    #sample_mask = np.zeros(num_samples, dtype='bool') # used samples

    logging.info(f"Sample Count {num_samples}")
    logging.info(f"Variant Count {num_vars}")

    # Okay, I accidentally implemented the greedy approach... 
    # So I'm curious how topN works, but until then, whatever
    out = open(out_fn, 'w')
    out.write("sample\tvarcount\tnew_count\ttot_captured\tpct_captured\n")
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
        new_variant_count = np.max(cur_sample_count)
        # don't want to use these variants anymore
        variant_mask = variant_mask | v_count[:, use_sample]
        # or this sample
        #sample_mask[use_sample] = True

        # our running total number of variants
        upto_now = variant_mask.sum()

        out.write(f"{vcf_samples[use_sample]}\t")
        out.write(f"{use_sample_variant_count}\t")
        out.write(f"{new_variant_count}\t")
        out.write(f"{upto_now}\t")
        out.write(f"{round(upto_now / num_vars, 4)}\n")

def load_files(in_files, safe=False):
    """
    Load and concatenate multiple files
    """
    logging.info(f"Loading {len(in_files)} files")
    samples = None
    parts = []
    for i in in_files:
        if i.endswith((".vcf.gz", ".vcf")):
            p = read_vcf(i)
        elif i.endswith(".jl"):
            p = joblib.load(i)
        if samples is None:
            samples = p['samples']
        elif safe:
            assert samples == p['samples']
        
        parts.append(p['GT'])
    logging.info("Concatenating")
    return {'GT':np.concatenate(parts), 'samples':samples}
            
def calc_main(cmdargs):
    args = parse_args(cmdargs)
    
    data = load_files(args.in_files, args.safe)
    calculate(data, args.out, args.keep)

    logging.info("Finished")

if __name__ == '__main__':
    calc_main(sys.argv[1])
