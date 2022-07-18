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
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args

def greedy_select(gt_matrix, vcf_samples, max_reporting, include, exclude, af, weights):
    """
    Greedy calculation
    yields list of samples
    """
    num_vars = gt_matrix.shape[0]
    # True where we've used the variant
    variant_mask = np.zeros(num_vars, dtype='bool')
    total_variant_count = gt_matrix.sum(axis=0)

    # get rid of exclude up front
    logging.info(f"Excluding {len(exclude)} samples")
    # False where we've used the sample
    sample_mask = ~np.isin(vcf_samples, exclude)

    tot_captured = 0

    logging.info(f"Including {len(include)} samples")
    for inc in include:
        use_sample = np.where(vcf_samples == inc)[0][0]
        variant_count = total_variant_count[use_sample]
        new_variant_count = gt_matrix[~variant_mask].sum(axis=0)[use_sample]
        variant_mask = variant_mask | gt_matrix[:, use_sample]
        sample_mask[use_sample] = False
        tot_captured  += new_variant_count
        yield [vcf_samples[use_sample], variant_count, new_variant_count,
               tot_captured, round(tot_captured / num_vars, 4)]

    for _ in range(max_reporting - len(include)):
        # of the variants remaining
        cur_view = gt_matrix[~variant_mask]
        # how many variants per sample
        cur_sample_count = cur_view.sum(axis=0) * sample_mask

        # scoring manipulations
        sample_scores = None
        if af is None and weights is None:
            sample_scores = cur_sample_count
        else:
            if af is not None:
                sample_scores = (cur_view * af[~variant_mask]).sum(axis=0) * sample_mask
            else:
                # Need to work on a copy so the values aren't destroyed
                sample_scores = cur_sample_count.copy()
            if weights is not None:
                sample_scores = sample_scores * weights

        # use highest score
        use_sample = np.argmax(sample_scores)
        # how many does this sample have overall?
        variant_count = total_variant_count[use_sample]
        # number of new variants added
        new_variant_count = cur_sample_count[use_sample]
        # don't want to use these variants anymore
        variant_mask = variant_mask | gt_matrix[:, use_sample]
        # or this sample
        sample_mask[use_sample] = False
        # update running total number of variants
        tot_captured += new_variant_count
        # stop running if we're out of new variants
        if new_variant_count == 0:
            logging.warning("Ran out of new variants")
            break

        yield [vcf_samples[use_sample], variant_count, new_variant_count,
               tot_captured, round(tot_captured / num_vars, 4)]

def run_selection(data, out_fn, max_reporting=0.02, include=None, exclude=None, weights=None):
    """
    Setup the selection calculation
    if max_reporting [0,1], select that percent of samples
    if max_reporting >= 1, select that number of samples
    """
    gt_matrix = data['GT']
    vcf_samples = data['samples']
    af_data = data["AF"]

    num_samples = gt_matrix.shape[1]
    num_vars = gt_matrix.shape[0]

    max_reporting = int(num_samples * max_reporting) if max_reporting < 1 else int(max_reporting)

    logging.info(f"Sample Count {num_samples}")
    logging.info(f"Variant Count {num_vars}")

    sample_weights = None
    if weights is not None:
        sample_weights = np.ones(num_samples)
        for pos, i in enumerate(vcf_samples):
            if i in weights.index:
                sample_weights[pos] = weights.loc[i]

    with open(out_fn, 'w') as out:
        out.write("sample\tvar_count\tnew_count\ttot_captured\tpct_captured\n")
        m_iter = greedy_select(gt_matrix, vcf_samples, max_reporting, include,
                               exclude, af_data, sample_weights)
        for result in m_iter:
            out.write("\t".join([str(_) for _ in result]) + '\n')

def samp_same(a, b):
    """
    Make sure samples are identical
    return true if they're identical
    """
    return len(a) == len(b) and np.equal(a, b).all()

def load_files(in_files, lowmem=False, load_af=False):
    """
    Load and concatenate multiple files
    """
    logging.info(f"Loading {len(in_files)} files")
    samples = None
    gt_parts = []
    af_parts = []
    for i in in_files:
        if i.endswith((".vcf.gz", ".vcf")):
            dat = read_vcf(i, lowmem, load_af)
        elif i.endswith(".jl"):
            dat = joblib.load(i)

        if samples is None:
            samples = dat['samples']
        elif not samp_same(samples, dat['samples']):
            logging.critical(f"Different samples in {i}")
            sys.exit(1)

        if dat['packedbits']:
            upack = np.unpackbits(dat['GT'], axis=1, count=len(dat['samples']))
            gt_parts.append(upack.astype(bool))
        else:
            gt_parts.append(dat['GT'])

        af_parts.append(dat['AF'])

    if len(gt_parts) > 1:
        logging.info("Concatenating")
        return  {'GT':np.concatenate(gt_parts),
                 'samples':samples,
                 'AF': np.concatenate(af_parts) if load_af else None}

    return  {'GT':gt_parts[0],
             'samples':samples,
             'AF': af_parts[0] if load_af else None}

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

    run_selection(data, args.out, args.count,
                  args.include, args.exclude,
                  args.weights)

    logging.info("Finished")
