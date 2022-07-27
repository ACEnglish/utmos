"""
Select fewest samples with maximum number of variants
"""
import os
import sys
import logging
import argparse
import tempfile

import h5py
import joblib
import truvari
import numpy as np
import pandas as pd

from utmos.convert import read_vcf
MAXMEM=2 # in GB

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

def do_lowmem_summation(matrix, variant_mask, sample_mask): #pylint:disable=unused-argument
    """
    Sum the matrix along axis=0
    if the matrix is an h5py.File, we'll operate in chunks
    otherwise, it is all in memory and simple to do
    """
    c_size = matrix.chunks
    logging.debug("csize is %s for shape of %s", c_size, matrix.shape)

    #I'm assuming chunks are full rows... may not be safe.
    n_rows, n_cols = matrix.shape

    one_row_size = n_cols * 4 / 1e9
    per_chunk_size = one_row_size * c_size[0]

    max_row_cnt = max(1, max(1, MAXMEM) / one_row_size)
    if max_row_cnt < c_size[0]:
        max_row_cnt = c_size[0]
    # round to nearest chunk boundary
    partition = max_row_cnt % c_size[0]
    if partition == 0:
        pass
    elif partition >= c_size[0] / 2:
        # round up
        max_row_cnt += c_size[0] - partition
    else:
        max_row_cnt -= partition
    max_row_cnt = min(int(max_row_cnt), n_rows)

    m_sum = np.zeros(n_cols)
    for row_start in range(0, n_rows, max_row_cnt):
        m_slice = slice(row_start, row_start + max_row_cnt, 1)
        logging.debug("slicing %d rows between %d:%d ~%.3fGB", max_row_cnt, m_slice.start, m_slice.stop, per_chunk_size)
        m_sum += matrix[m_slice, :].sum(axis=0)
    m_sum *= sample_mask
    return m_sum


def calculate_scores(gt_matrix, variant_mask, sample_mask, af_matrix, sample_weights, summer=do_summation):
    """
    return number of variants and the scores per-sample for this iteration
    summer is the method to do the summation (experimental)
    """
    logging.debug("gt_matrix sum")
    cur_sample_count = summer(gt_matrix, variant_mask, sample_mask)

    sample_scores = cur_sample_count
    if af_matrix is not None:
        logging.debug("af_matrix sum")
        sample_scores = summer(af_matrix, variant_mask, sample_mask)
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
        total_variant_count = do_summation(gt_matrix, variant_mask, sample_mask)
    else:
        total_variant_count = gt_matrix.sum(axis=0)

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

        yield [vcf_samples[use_sample], int(variant_count), int(new_variant_count),
               int(tot_captured), round(tot_captured / num_vars, 4)]

def rewrite_smaller_hdf5(gt_matrix, af_matrix, variant_mask, sample_mask, temp_name):
    """
    Attempting to update gt_matrix and af_matrix by pulling all the data into a smaller dataframe
    """
    # should remove previous temp
    logging.debug('rewriting to %s', temp_name)
    n_cols = sample_mask.sum()
    n_rows = (~variant_mask).sum()
    c_rows = min(n_rows, max(1, int(1e6 / 4 / n_cols)))
    c_size = (c_rows, n_cols)
    if 0 in c_size:
        logging.debug("finished? %s", c_size)
        return None, None
    with h5py.File(temp_name, 'w') as hf:
        dset = hf.create_dataset('GT', shape=(n_rows, n_cols), compression="lzf", dtype='bool', chunks=c_size)
        m_tmp = gt_matrix[~variant_mask, :]
        m_tmp = np.ascontiguousarray(m_tmp[:, sample_mask][:])
        dset.write_direct(m_tmp)
        if af_matrix:
            dset = hf.create_dataset('AF_matrix', shape=(n_rows, n_cols), compression="lzf", dtype='float', chunks=c_size)
            m_tmp = af_matrix[~variant_mask, :]
            m_tmp = np.ascontiguousarray(m_tmp[:, sample_mask][:])
            dset.write_direct(m_tmp)

    logging.debug('reloading')
    new_fh = h5py.File(temp_name, 'r')
    gt_matrix = new_fh["GT"]
    af_matrix = new_fh["AF_matrix"] if af_matrix else None
    return gt_matrix, af_matrix

def under_max_mem_usage(shape, with_af=False):
    """
    Simple memory usage estimate in GB
    """
    return (shape[0] * shape[1] * 4 * (int(with_af) + 1)) / 1e9 < MAXMEM

# pylint:disable=too-many-locals
def greedy_mem_select(gt_matrix, select_count, vcf_samples, variant_mask, sample_mask, af_matrix=None,
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

    Expects input matrices to be h5py Datasets. Will Do an iterative rewrite
    """
    if not isinstance(gt_matrix, h5py.Dataset):
        logging.error("Can only run greedy_mem on h5df")
        sys.exit(1)
    logging.debug("running greedy_mem mode")

    num_vars = gt_matrix.shape[0]

    # Only need to calculate this once
    logging.debug("getting total_variant_count")
    if isinstance(gt_matrix, h5py.Dataset):
        total_variant_count = do_lowmem_summation(gt_matrix, variant_mask, sample_mask)
    else:
        total_variant_count = gt_matrix.sum(axis=0)

    logging.debug("per-sample variant mean %.1f", total_variant_count.mean())

    tot_captured = 0
    prev_tmp_name = None
    for _ in range(select_count):
        cur_sample_count, sample_scores = calculate_scores(gt_matrix, variant_mask, sample_mask, af_matrix,
                                                           sample_weights, do_lowmem_summation)

        # use highest score
        use_sample = np.argmax(sample_scores)
        use_sample_name = vcf_samples[use_sample]
        # how many does this sample have overall?
        variant_count = total_variant_count[use_sample]
        # number of new variants added
        new_variant_count = cur_sample_count[use_sample]
        # don't want to use these variants anymore - this is anti-pattern to chunks
        sample_mask[use_sample] = False
        # update running total number of variants
        tot_captured += new_variant_count
        # stop running if we're out of new variants
        if new_variant_count == 0:
            logging.warning("Ran out of new variants")
            return

        # Gonna be dropping this sample
        vcf_samples = vcf_samples[sample_mask]
        total_variant_count = total_variant_count[sample_mask]
        if sample_weights is not None:
            sample_weights = sample_weights[sample_mask]

        variant_mask = gt_matrix[:, use_sample].astype(bool)

        ## Re-writing
        pre_shape = gt_matrix.shape
        temp_name = os.path.join(tempfile.gettempdir(),
                                 next(tempfile._get_candidate_names()) + '.utmos.tmp.hdf5') # pylint:disable=protected-access, stop-iteration-return

        gt_matrix, af_matrix = rewrite_smaller_hdf5(gt_matrix, af_matrix, variant_mask, sample_mask, temp_name)

        # clean tmp files
        if prev_tmp_name is not None:
            logging.debug("removing %s", prev_tmp_name)
            os.remove(prev_tmp_name)
        prev_tmp_name = temp_name

        yield [use_sample_name, int(variant_count), int(new_variant_count),
               int(tot_captured), round(tot_captured / num_vars, 4)]
        if not gt_matrix: # Ran out of data
            return

        logging.debug("shape %s -> %s", pre_shape, gt_matrix.shape)

        # Hacky
        sample_mask = np.ones(gt_matrix.shape[1], dtype='bool')

        # Stop doing hdf5 rewrites when we can
        if isinstance(gt_matrix, h5py.Dataset) and under_max_mem_usage(gt_matrix.shape, af_matrix is not None):
            logging.info("Dataset small enough to hold in memory")
            gt_matrix = gt_matrix[:]
            af_matrix = af_matrix[:] if af_matrix else af_matrix
            for i in greedy_select(gt_matrix, select_count, vcf_samples, variant_mask, sample_mask, af_matrix,
                                   sample_weights):
                yield i
            return


SELECTORS = {"greedy": greedy_select,
             "greedy_mem": greedy_mem_select}
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

    gt_matrix = data['GT']
    af_matrix = None
    if af and isinstance(data, h5py.File):
        af_matrix = data["AF_matrix"]
    elif af:
        logging.info("Calculating AF matrix")
        af_matrix = gt_matrix * data['AF']
        af_matrix[np.isnan(af_matrix)] = 0

    #if logging.root.level <= logging.DEBUG:
    #    logging.debug("pre-flight variant check")
    #    cnt = gt_matrix[:,sample_mask].any(axis=1).sum()
    #    logging.debug("have %d variants to process", cnt)

    m_select = SELECTORS[mode]
    if isinstance(data, h5py.File) and under_max_mem_usage(gt_matrix.shape, af_matrix is not None):
        logging.info("Dataset small enough to hold in memory")
        gt_matrix = gt_matrix[:]
        af_matrix = af_matrix[:] if af_matrix else af_matrix
        return greedy_select(gt_matrix, select_count, vcf_samples, variant_mask, sample_mask, af_matrix, sample_weights)

    return m_select(gt_matrix, select_count, vcf_samples, variant_mask, sample_mask, af_matrix,
                    sample_weights)

def write_append_hdf5(cur_part, out_name, is_first=False):
    """
    Handle the hdf5 file
    """
    # Future - make af_matrix optional
    logging.debug('calc af')

    #logging.debug("sorting by AF") (doesn't make sense with greedy_mem not slicing as much anymore)
    #m_order = cur_part["AF"].argsort(axis=0)
    #cur_part["GT"] = cur_part["GT"][m_order]
    #cur_part["AF"] = cur_part["AF"][m_order]

    af_matrix = cur_part["GT"] * cur_part["AF"]
    af_matrix[np.isnan(af_matrix)] = 0
    n_cols = cur_part['GT'].shape[1]
    c_size = (max(1, int(1e6 / 4 / n_cols)), n_cols)
    if is_first:
        logging.debug('write')
        with h5py.File(out_name, 'w') as hf:
            hf.create_dataset('GT', data=cur_part["GT"], compression="lzf", chunks=c_size, maxshape=(None, n_cols))
            hf.create_dataset('AF_matrix', data=af_matrix, compression="lzf", chunks=c_size, maxshape=(None, n_cols))
            hf.create_dataset('samples', data=cur_part["samples"], compression="lzf", chunks=True, maxshape=(None,))
        return

    with h5py.File(out_name, 'a') as hf:
        logging.debug('append')
        hf["GT"].resize((hf["GT"].shape[0] + cur_part["GT"].shape[0]), axis = 0)
        hf["GT"][-cur_part["GT"].shape[0]:] = cur_part["GT"]

        hf["AF_matrix"].resize((hf["AF_matrix"].shape[0] + af_matrix.shape[0]), axis = 0)
        hf["AF_matrix"][-af_matrix.shape[0]:] = af_matrix


#pylint:disable=too-many-statements
def load_files(in_files, lowmem=None, buffer=32768):
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

        upack = np.unpackbits(dat['GT'], axis=1, count=len(dat['samples'])).astype(bool)
        uninf_filter = upack.any(axis=1)
        logging.debug("fitering %d uninformative variants", (~uninf_filter).sum())

        gt_parts.append(upack[uninf_filter])
        af_parts.append(dat['AF'][uninf_filter])
        m_count = dat["GT"].shape[0]

        load_count += 1
        load_row_count += m_count
        load_buffer_count += m_count
        if lowmem is not None and (load_buffer_count >= buffer or pos == len(in_files) - 1):
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
            # reset buffering
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
    parser.add_argument("-c", "--count", type=float, default=0.02,
                        help="Number of samples to select as a percent if <1 or count if >=1 (%(default)s)")
    parser.add_argument("-o", "--out", type=str, default="/dev/stdout",
                        help="Output file (stdout)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")

    scoreg = parser.add_argument_group("Scoring Arguments")
    scoreg.add_argument("--af", action="store_true",
                        help="Weigh variants by allele frequency")
    scoreg.add_argument("--weights", type=str, default=None,
                        help="Tab-delimited file of sample weights")
    scoreg.add_argument("--subset", type=str, default=None,
                        help="Filename with or Comma-separated list of samples to analyze")
    scoreg.add_argument("--exclude", type=str, default=None,
                        help="Filename with or Comma-separated list of samples to exclude selection")

    mperfg = parser.add_argument_group("Memory Arguments")
    mperfg.add_argument("--lowmem", type=str, default=None,
                        help="Name of concatenated hdf5 file to create/use (%(default)s)")
    mperfg.add_argument("--buffer", type=int, default=32768,
                        help="When using `--lowmem`, number of variants to buffer during concatenation (%(default)s)")
    mperfg.add_argument("--maxmem", type=int, default=2,
                        help="Maximum amount of memory in (GB) to use with --lowmem. 0 keeps data in hdf5 (%(default)s)")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


def select_main(cmdargs):
    """
    Main
    """
    global MAXMEM
    args = parse_args(cmdargs)

    data = load_files(args.in_files, args.lowmem, args.buffer)
    args.subset = parse_sample_lists(args.subset)
    #args.include = parse_sample_lists(args.include)
    args.exclude = parse_sample_lists(args.exclude)
    args.weights = parse_weights(args.weights)
    mode = 'greedy_mem' if args.lowmem else 'greedy'
    MAXMEM = args.maxmem
    with open(args.out, 'w') as fout:
        fout.write("sample\tvar_count\tnew_count\ttot_captured\tpct_captured\n")
        m_iter = run_selection(data, args.count, mode, args.subset,
                               args.exclude, args.af, args.weights)
        for result in m_iter:
            logging.info("Selected %s (%s)", result[0], result[4])
            fout.write("\t".join([str(_) for _ in result]) + '\n')

    logging.info("Finished")
