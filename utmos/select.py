"""
Select fewest samples with maximum number of variants
"""
import os
import sys
import logging
import argparse

import h5py
import joblib
import truvari
import numpy as np
import pandas as pd

from utmos.convert import read_vcf

def do_summation(matrix, variant_mask, sample_mask):
    """
    Sum the matrix along axis=0
    if the matrix is an h5py.File, we'll operate in chunks
    otherwise, it is all in memory and simple to do
    """
    if not isinstance(matrix, h5py.Dataset):
        return matrix[~variant_mask].sum(axis=0) * sample_mask

    m_sum = np.zeros(matrix.shape[1])
    for i in matrix.iter_chunks():
        cur_v_mask = variant_mask[i[0]]
        cur_chunk = matrix[i]
        m_sum[i[1]] += cur_chunk[~cur_v_mask].sum(axis=0) * sample_mask[i[1]]
    return m_sum

def calculate_scores(gt_matrix, variant_mask, sample_mask, af_matrix, sample_weights):
    """
    return number of variants and the scores per-sample for this iteration
    """
    logging.debug("gt_matrix sum")
    cur_sample_count = do_summation(gt_matrix, variant_mask, sample_mask)

    sample_scores = cur_sample_count
    if af_matrix is not None:
        logging.debug("af_matrix sum")
        sample_scores = do_summation(af_matrix, variant_mask, sample_mask)
    elif sample_weights is not None:
        sample_scores = cur_sample_count.copy()
    if sample_weights is not None:
        logging.debug("applying weights")
        sample_scores = sample_scores * sample_weights
    return cur_sample_count, sample_scores

def greedy_select(gt_matrix, select_count, vcf_samples, variant_mask, sample_mask, af_matrix=None,
                  sample_weights=None):
    """
    Greedy calculation
    yields rows of each selected sample's information

    gt_matrix = boolean matrix of genotype presence
    select_count = how many samples we'll be selecting
    vcf_samples = identifiers of sample names (len == gt_matrix.shape[1])
    variant_mask = boolean matrix of variants where True == used
    sample_mask = boolean matrix of samples where True == use
    af_matrix = (optional) the af_matrix scores (af_matrix.shape == gt_matrix.shape)
    sample_weights = (optional) the weights to apply to each iteration's sample.sum (len == gt_matrix.shape[0])
    """
    num_vars = gt_matrix.shape[0]
    # Only need to calculate this once
    logging.debug("getting total_variant_count")
    if isinstance(gt_matrix, h5py.Dataset):
        v_mask = np.zeros(gt_matrix.shape[0], dtype='bool')
        s_mask = np.ones(gt_matrix.shape[1], dtype='bool')
        total_variant_count = do_summation(gt_matrix, v_mask, s_mask)
    else:
        total_variant_count = gt_matrix.sum(axis=0)

    logging.debug("per-sample variant mean %.1f", total_variant_count.mean())

    tot_captured = 0
    for _ in range(select_count):
        cur_sample_count, sample_scores = calculate_scores(gt_matrix, variant_mask, sample_mask, af_matrix,
                                                           sample_weights)

        # use highest score
        use_sample = np.argmax(sample_scores)
        # how many does this sample have overall?
        variant_count = total_variant_count[use_sample]
        # number of new variants added
        new_variant_count = cur_sample_count[use_sample]
        # don't want to use these variants anymore - this is anti-pattern to chunks
        logging.debug('updating variant mask') 
        variant_mask = variant_mask | gt_matrix[:, use_sample]
        # or this sample
        sample_mask[use_sample] = False
        # update running total number of variants
        tot_captured += new_variant_count
        # stop running if we're out of new variants
        if new_variant_count == 0:
            logging.warning("Ran out of new variants")
            break

        yield [vcf_samples[use_sample], int(variant_count), int(new_variant_count),
               int(tot_captured), round(tot_captured / num_vars, 4)]

def topN_select(gt_matrix, vcf_samples, max_reporting, include, exclude, af, weights, is_random=False):
    """
    Select the topN samples
    I reuse this method to do random selection as well
    """
    num_vars = gt_matrix.shape[0]
    # True where we've used the variant
    variant_mask = np.zeros(num_vars, dtype='bool')

    sample_frame = pd.DataFrame(pd.Series(gt_matrix.sum(axis=0)), columns=['count'])
    sample_frame["sample"] = vcf_samples

    # scoring
    sort_key = 'count'
    if af is not None:
        sample_frame['score'] = (gt_matrix * af).sum(axis=0)
        sort_key = 'score'
    elif weights is not None:
        # Need to let weights work on a copy of counts so values aren't destroyed
        sample_frame['score'] = sample_frame["count"]
    if weights is not None:
        sort_key = 'score'
        sample_frame['score'] = sample_frame['score'] * weights

    logging.info(f"Excluding {len(exclude)} samples")
    logging.info(f"Including {len(include)} samples")

    tot_captured = 0
    for inc in include:
        exclude.append(inc)
        view = sample_frame[sample_frame['sample'] == inc].iloc[0]
        new_variant_count = gt_matrix[~variant_mask, view.name].sum(axis=0)
        tot_captured += new_variant_count
        variant_mask = variant_mask | gt_matrix[:, view.name]
        yield [inc, view['count'], new_variant_count,
               tot_captured, round(tot_captured / num_vars, 4)]

    max_reporting -= len(include)
    sample_frame.sort_values(sort_key, ascending=False, inplace=True)
    sample_mask = sample_frame[~sample_frame["sample"].isin(exclude)].index
    if is_random:
        m_iter = sample_frame.loc[sample_mask].sample(max_reporting)
    else:
        m_iter = sample_frame.loc[sample_mask][:max_reporting]

    for idx, d in m_iter.iterrows():
        new_variant_count = gt_matrix[~variant_mask, idx].sum()
        tot_captured += new_variant_count
        variant_mask = variant_mask | gt_matrix[:, idx]
        yield [d["sample"], d["count"], new_variant_count,
               tot_captured, round(tot_captured / num_vars, 4)]

def random_select(*args, **kwargs):
    """
    Call topN select in random mode
    """
    return topN_select(*args, **kwargs, is_random=True)

SELECTORS = {"greedy": greedy_select}
             # deprecated for now
             #"topN": topN_select,
             #"random": random_select}

def run_selection(data, select_count=0.02, mode='greedy', subset=None, exclude=None, af=False,
                  weights=None):
    """
    Setup the selection calculation
    if select_count [0,1], select that percent of samples
    if select_count >= 1, select that number of samples
    returns the generator created by the specified mode
    """
    num_vars, num_samples = data["GT"].shape
    logging.info("Sample Count %d", num_samples)
    logging.info("Variant Count %d", num_vars)

    select_count = max(1, int(num_samples * select_count) if select_count < 1 else int(select_count))
    logging.info("Selecting %d", select_count)

    # Build masks
    variant_mask = np.zeros(num_vars, dtype='bool')

    vcf_samples = data['samples']
    if isinstance(data, h5py.File):
        vcf_samples = data['samples'][:].astype(str)
    else:
        vcf_samples = vcf_samples.astype(str)

    sample_mask = np.ones(num_samples, dtype='bool')
    if subset:
        sample_mask = np.isin(vcf_samples, subset)
        logging.info("Subsetting to %d of %d samples", sample_mask.sum(), num_samples)
    if exclude:
        logging.info(f"Excluding {len(exclude)} samples")
        sample_mask = sample_mask & ~np.isin(vcf_samples, exclude)
    logging.debug(f"ending with {sample_mask.sum()} samples")
       
    sample_weights = None
    if weights is not None:
        logging.info("Setting %d weights", len(weights))
        sample_weights = np.ones(num_samples)
        for pos, i in enumerate(vcf_samples):
            if i in weights.index:
                sample_weights[pos] = weights.loc[i]
    logging.debug("finished with weights")

    gt_matrix = data['GT']
    af_matrix = None
    if af and isinstance(data, h5py.File):
        af_matrix = data["AF_matrix"]
    elif af:
        logging.info("Calculating AF matrix")
        af_matrix = gt_matrix * data['AF']

    #if logging.root.level <= logging.DEBUG:
    #    logging.debug("pre-flight variant check")
    #    cnt = gt_matrix[:,sample_mask].any(axis=1).sum()
    #    logging.debug("have %d variants to process", cnt)

    m_select = SELECTORS[mode]
    return m_select(gt_matrix, select_count, vcf_samples, variant_mask, sample_mask, af_matrix,
                    sample_weights)

def samp_same(a, b):
    """
    Make sure samples are identical
    return true if they're identical
    """
    return len(a) == len(b) and np.equal(a, b).all()

def write_append_hdf5(cur_part, out_name, is_first=False):
    """
    Handle the hdf5 file
    """
    # Future - make af_matrix optional
    logging.debug('calc af')
    af_matrix = cur_part["GT"] * cur_part["AF"]
    n_cols = cur_part['GT'].shape[1]
    c_size = (max(1, int(1e6 / 4 / n_cols)), n_cols)
    if is_first:
        logging.debug('write')
        with h5py.File(out_name, 'w') as hf:
            hf.create_dataset('GT', data=cur_part["GT"], compression="lzf", chunks=c_size, maxshape=(None, n_cols))
            #hf.create_dataset('AF', data=cur_part["AF"], compression="lzf", chunks=True, maxshape=(None, 1))
            hf.create_dataset('AF_matrix', data=af_matrix, compression="lzf", chunks=c_size, maxshape=(None, n_cols))
            hf.create_dataset('samples', data=cur_part["samples"], compression="lzf", chunks=True, maxshape=(None,))
        return

    with h5py.File(out_name, 'a') as hf:
        logging.debug('append')
        hf["GT"].resize((hf["GT"].shape[0] + cur_part["GT"].shape[0]), axis = 0)
        hf["GT"][-cur_part["GT"].shape[0]:] = cur_part["GT"]

        #hf["AF"].resize((hf["AF"].shape[0] + cur_part["AF"].shape[0]), axis = 0)
        #hf["AF"][-cur_part["AF"].shape[0]:] = cur_part["AF"]

        hf["AF_matrix"].resize((hf["AF_matrix"].shape[0] + af_matrix.shape[0]), axis = 0)
        hf["AF_matrix"][-af_matrix.shape[0]:] = af_matrix


def load_files(in_files, lowmem=None, chunk_length=32768):
    """
    Load and concatenate multiple files
    if lowmem is provided, in_files are concatenated into an hdf5. This may be a little slower,
    but will use less memory since it won't need to hold two copies of the
    data in memory during the concatenation step.
    If the af_matrix is going to be used, lowmem will also create that per-in_file and write it to the hdf5 file
    """
    file_cnt = len(in_files)
    logging.info(f"Loading {file_cnt} files")
    samples = None
    gt_parts = []
    af_parts = []
    load_count = 0
    load_row_count = 0
    load_buffer_count = 0
    is_first = True
    for pos, i in enumerate(in_files):
        if i.endswith((".vcf.gz", ".vcf")):
            dat = read_vcf(i, lowmem is not None)
        elif i.endswith(".jl"):
            dat = joblib.load(i)
        elif i.endswith(".hdf5"):
            if len(in_files) > 1:
                logging.error("Only one hdf5 file at a time!")
                sys.exit(1)
            logging.info("Loading single hdf5")
            lowmem = i
            break
        else:
            logging.error("Unknown filetype %s. Expected `.vcf[.gz]`, `.jl`, or `.hdf5`", i)
            sys.exit(1)

        if samples is None:
            samples = dat['samples'].astype('S')
        elif not samp_same(samples, dat['samples'].astype('S')):
            logging.critical(f"Different samples in {i}")
            sys.exit(1)

        upack = np.unpackbits(dat['GT'], axis=1, count=len(dat['samples']))
        gt_parts.append(upack.astype(bool))
        af_parts.append(dat['AF'])
        m_count = dat["GT"].shape[0]

        load_count += 1
        load_row_count += m_count
        load_buffer_count += m_count
        if lowmem is not None and (load_buffer_count >= chunk_length or pos == len(in_files) - 1):
            logging.debug("dumping chunk with %d vars", load_buffer_count)
            cur_chunk = None
            if len(gt_parts) > 1:
                cur_chunk = {'GT': np.concatenate(gt_parts),
                             'samples': samples,
                             'AF': np.concatenate(af_parts)}
            else:
                cur_chunk = {'GT':gt_parts[0],
                             'samples':samples,
                             'AF': af_parts[0]}
            
            write_append_hdf5(cur_chunk, lowmem, is_first)
            # reset
            load_buffer_count = 0
            is_first = False
            gt_parts = []
            af_parts = []

        logging.debug("Loaded %d of %d (%.2f%%) with %d vars (%d buff)", load_count, file_cnt,
                      load_count / file_cnt * 100, load_row_count, load_buffer_count)

    if lowmem is not None:
        return h5py.File(lowmem, 'r')

    if len(gt_parts) > 1:
        logging.info("Concatenating")
        return  {'GT':np.concatenate(gt_parts),
                 'samples':samples,
                 'AF': np.concatenate(af_parts)}

    return  {'GT':gt_parts[0],
             'samples':samples,
             'AF': af_parts[0]}

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

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="select", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_files", nargs="+", type=str,
                        help="Input VCF or jl files")
    parser.add_argument("--lowmem", type=str, default=None,
                        help="Name of concatenated hdf5 file to create, which reduces memory usage (%(default)s)")
    parser.add_argument("--chunk-length", type=int, default=32768,
                        help="When using `--lowmem`, number of variants to process at a time (%(default)s)")
    parser.add_argument("-o", "--out", type=str, default="/dev/stdout",
                        help="Output file (stdout)")
    parser.add_argument("-m", "--mode", type=str, default='greedy', choices=SELECTORS.keys(),
                        help="Selection algo to use greedy (default), topN, or random")
    parser.add_argument("-c", "--count", type=float, default=0.02,
                        help="Number of samples to select as a percent if <1 or count if >=1 (%(default)s)")
    parser.add_argument("--af", action="store_true",
                        help="Weigh variants by allele frequency")
    parser.add_argument("--weights", type=str, default=None,
                        help="Tab-delimited file of sample weights")
    parser.add_argument("--subset", type=str, default=None,
                        help="Filename with or Comma-separated list of sample subset to analyze")
    #parser.add_argument("--include", type=str, default=None,
                        #help="Filename with or Comma-separated list of samples to force selection")
    parser.add_argument("--exclude", type=str, default=None,
                        help="Filename with or Comma-separated list of samples to exclude selection")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


def select_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)

    data = load_files(args.in_files, args.lowmem, args.chunk_length)
    args.subset = parse_sample_lists(args.subset)
    #args.include = parse_sample_lists(args.include)
    args.exclude = parse_sample_lists(args.exclude)
    args.weights = parse_weights(args.weights)

    with open(args.out, 'w') as fout:
        fout.write("sample\tvar_count\tnew_count\ttot_captured\tpct_captured\n")
        m_iter = run_selection(data, args.count, args.mode, args.subset,
                               args.exclude, args.af, args.weights, args.chunk_length)
        for result in m_iter:
            logging.info("Selected %s (%s)", result[0], result[4])
            fout.write("\t".join([str(_) for _ in result]) + '\n')

    logging.info("Finished")
