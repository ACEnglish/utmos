"""
Select fewest samples with maximum number of variants
"""
import os
import sys
import logging
import argparse

import joblib
import truvari
import numpy as np
import pandas as pd

from utmos.convert import read_vcf

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="select", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_files", nargs="+", type=str,
                        help="Input VCF or jl files")
    parser.add_argument("--lowmem", action="store_true",
                        help="Use temporary files to lower memory usage during vcf conversion (%(default)s)")
    parser.add_argument("-o", "--out", type=str, default="/dev/stdout",
                        help="Output file (stdout)")
    parser.add_argument("-c", "--count", type=float, default=0.02,
                        help="Number of samples to select as a percent if <1 or count if >=1 (%(default)s)")
    parser.add_argument("--af", action="store_true",
                        help="Weigh variants by allele frequency")
    parser.add_argument("--weights", type=str, default=None,
                        help="Tab-delimited file of sample weights")
    parser.add_argument("--include", type=str, default=None,
                        help="Filename with or Comma-separated list of samples to force selection")
    parser.add_argument("--exclude", type=str, default=None,
                        help="Filename with or Comma-separated list of samples to exclude selection")
    args = parser.parse_args(args)
    truvari.setup_logging()
    return args

def greedy_calc(v_count, vcf_samples, max_reporting, include, exclude, af, weights):
    """
    Greedy calculation
    yields list of samples
    """
    num_vars = v_count.shape[0]
    variant_mask = np.zeros(num_vars, dtype='bool')
    # get rid of exclude up front
    logging.info(f"Excluding {len(exclude)} samples")
    sample_mask = ~np.isin(vcf_samples, exclude)

    logging.info(f"Including {len(include)} samples")
    # calculate this once
    total_variant_count = v_count.sum(axis=0)

    upto_now = 0
    for inc in include:
        use_sample = np.where(vcf_samples == inc)[0][0]
        cur_view = v_count[~variant_mask]
        cur_sample_count = cur_view.sum(axis=0)
        use_sample_variant_count = total_variant_count[use_sample]
        new_variant_count = cur_sample_count[use_sample]
        variant_mask = variant_mask | v_count[:, use_sample]
        upto_now += new_variant_count
        yield [vcf_samples[use_sample], use_sample_variant_count, new_variant_count,
               upto_now, round(upto_now / num_vars, 4)]

    for _ in range(max_reporting - len(include)): # picking the best N samples
        # of the variants remaining
        cur_view = v_count[~variant_mask]
        # how many variants per sample
        cur_sample_count = cur_view.sum(axis=0) * sample_mask

        # use the sample with the most variants
        # incorporate weights if needed
        cur_sample_weighted = None
        if af is not None:
            cur_sample_weighted = (cur_view * af[~variant_mask]).sum(axis=0)
        if weights is not None:
            if cur_sample_weighted is None:
                cur_sample_weighted = cur_sample_count.copy()
            cur_sample_weighted = cur_sample_weighted * weights

        if af is not None or weights is not None:
            use_sample = np.argmax(cur_sample_weighted)
        else:
            use_sample = np.argmax(cur_sample_count)

        # how many does this sample have overall?
        use_sample_variant_count = total_variant_count[use_sample]
        # number of new variants added
        new_variant_count = cur_sample_count[use_sample]
        # don't want to use these variants anymore
        variant_mask = variant_mask | v_count[:, use_sample]
        # or this sample
        sample_mask[use_sample] = False

        # our running total number of variants
        upto_now += new_variant_count
        # stop running if we're out of new variants
        if new_variant_count == 0:
            logging.warning("Ran out of new variants")
            break

        yield [vcf_samples[use_sample], use_sample_variant_count, new_variant_count,
               upto_now, round(upto_now / num_vars, 4)]

def calculate(data, out_fn, max_reporting=0.02, include=None, exclude=None, af=False, weights=None):
    """
    Do the selection calculation
    if max_reporting [0,1], select that percent of samples
    if max_reporting >= 1, select that number of samples
    """
    v_count = data['GT']
    vcf_samples = data['samples']
    af_data = None
    if af:
        af_data = data["AF"]
        af_data = af_data.reshape(af_data.shape[0], 1)

    num_samples = v_count.shape[1]
    num_vars = v_count.shape[0]

    max_reporting = int(num_samples * max_reporting) if max_reporting < 1 else int(max_reporting)

    logging.info(f"Sample Count {num_samples}")
    logging.info(f"Variant Count {num_vars}")

    sample_weights = None
    if weights is not None:
        sample_weights = np.zeros(num_samples) + 1
        for pos, i in enumerate(vcf_samples):
            if i in weights.index:
                sample_weights[pos] = weights.loc[i]

    with open(out_fn, 'w') as out:
        out.write("sample\tvar_count\tnew_count\ttot_captured\tpct_captured\n")
        for result in greedy_calc(v_count, vcf_samples, max_reporting, include, exclude, af_data, sample_weights):
            out.write("\t".join([str(_) for _ in result]) + '\n')

def samp_same(a, b):
    """
    Make sure samples are identical
    return true if they're identical
    """
    return len(a) == len(b) and np.equal(a, b).all()

def load_files(in_files, lowmem=False, af=False):
    """
    Load and concatenate multiple files
    """
    logging.info(f"Loading {len(in_files)} files")
    samples = None
    gt_parts = []
    af_parts = []
    for i in in_files:
        if i.endswith((".vcf.gz", ".vcf")):
            p = read_vcf(i, lowmem, af)
        elif i.endswith(".jl"):
            p = joblib.load(i)
        if samples is None:
            samples = p['samples']
        elif not samp_same(samples, p['samples']):
            logging.critical(f"Different samples in {i}")
            sys.exit(1)

        gt_parts.append(p['GT'])
        if af:
            af_parts.append(p['AF'])

    ret = None
    if len(gt_parts) > 1:
        logging.info("Concatenating")
        ret =  {'GT':np.concatenate(gt_parts),
                'samples':samples}
        if af:
            ret['AF'] = np.concatenate(af_parts)
    else:
        ret =  {'GT':gt_parts[0],
                'samples':samples}
        if af:
            ret['AF'] = af_parts
        
    return ret

def parse_sample_lists(argument):
    """
    Parse the --include or --exclude argument.
    """
    ret = []
    if not argument:
        return ret
    if os.path.exists(argument):
        with open(argument, 'r') as fh:
            ret = [_.strip() for _ in fh]
    else:
        ret = argument.split(",")
    return ret

def parse_weights(argument):
    """
    Parse the weights file
    """
    if not argument:
        return None
    data = pd.read_csv(argument, sep='\t', header=None)
    data.columns = ['sample', 'weight']
    data.set_index('sample', inplace=True)
    return data

def select_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)

    data = load_files(args.in_files, args.lowmem, args.af)
    args.include = parse_sample_lists(args.include)
    args.exclude = parse_sample_lists(args.exclude)
    args.weights = parse_weights(args.weights)

    calculate(data, args.out, args.count,
              args.include, args.exclude,
              args.af, args.weights)

    logging.info("Finished")
