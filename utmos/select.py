"""
Select fewest samples with maximum number of variants
"""
import os
import sys
import json
import logging
import argparse

import h5py
import joblib
import truvari
import numpy as np
import pandas as pd

from utmos.convert import read_vcf

MAXMEM = 2  # in GB


#############
# Core code #
#############
def do_sum(matrix, sample_mask):
    m_sum = np.zeros(matrix.shape[1])
    m_tot = np.zeros(matrix.shape[1])
    c_mask = np.where(~sample_mask)
    for row in matrix:
        if row[c_mask].any():
            continue
        m_sum += row
        if matrix.dtype == bool:
            m_tot += row
        else:
            m_tot += row != 0
    m_sum *= sample_mask
    m_tot *= sample_mask
    return m_sum, m_tot


def calculate_scores(matrix, sample_mask, sample_weights):
    """
    calculate the best scoring sample,
    sumfunc is the method to do matrix summation

    updates sample_mask in place
    returns tuple of:
        column index of the highest score
        new_row_count for highest score column index
    """
    sample_scores, cur_sample_count = do_sum(matrix, sample_mask)

    if sample_weights is not None:
        logging.debug("applying weights")
        sample_scores = sample_scores * sample_weights

    use_sample = np.argmax(sample_scores)
    new_variant_count = cur_sample_count[use_sample]

    sample_mask[use_sample] = False

    return use_sample, new_variant_count


def is_memsafe(shape, with_af=False):
    """
    Simple memory usage estimate in GB
    returns true if we think we can hold a full dataset in memory
    """
    data_size = (shape[0] * shape[1] * 4 * (int(with_af) + 1)) / 1e9
    logging.debug("Estimated data size %.2fGB", data_size)
    return data_size < MAXMEM


##############
# Algorithms #
##############
def greedy_select(matrix,
                  total_variant_count,
                  select_count,
                  vcf_samples,
                  sample_mask,
                  sample_weights=None):
    """
    Greedy calculation
    yields rows of each selected sample's information

    matrix:              genotype data
    total_variant_count: total number of variants per-sample
    select_count:        how many samples we'll be selecting
    variant_mask:        boolean matrix of variants where True == used
    sample_mask:         boolean matrix of samples where True == use
    sample_weights:      (optional) the weights to apply to each iteration's sample.sum (len == gt_matrix.shape[0])

    Expects input matrices to be h5py Datasets.
    Will do an iterative rewrite until is_memsafe == True
    """
    num_vars = matrix.shape[0]

    tot_captured = 0
    for _ in range(select_count):
        use_sample, new_variant_count = calculate_scores(matrix, sample_mask, sample_weights)

        use_sample_name = vcf_samples[use_sample]
        variant_count = total_variant_count[use_sample]
        tot_captured += new_variant_count

        yield [
            use_sample_name,
            int(variant_count),
            int(new_variant_count),
            int(tot_captured),
            round(tot_captured / num_vars, 4)
        ]

        if tot_captured >= num_vars:
            logging.warning("Ran out of new variants")
            return

        # can mem? put it in
        # need to change shape by how many we could mask
        if isinstance(matrix, h5py.Dataset):
            n_var = tot_captured
            n_samp = sample_mask.sum()
            if is_memsafe((n_var, n_samp)):
                logging.info("Dataset small enough to hold in memory")
                n_matrix = np.zeros((n_var, n_samp), dtype=matrix.dtype)
                m_pos = 0
                c_mask = np.where(~sample_mask)
                for row in matrix:
                    if row[c_mask].any():
                        continue
                    n_matrix[m_pos] = row[sample_mask]
                    m_pos += 1
                # Drop used samples
                vcf_samples = vcf_samples[sample_mask]
                total_variant_count = total_variant_count[sample_mask]
                if sample_weights is not None:
                    sample_weights = sample_weights[sample_mask]
                sample_mask = np.ones(n_matrix.shape[1], dtype='bool')
                matrix = n_matrix


# deprecated for now
#"topN": topN_select,
#"random": random_select}

####################
# Setup/Management #
####################
def run_selection(data, select_count=0.02, subset=None, exclude=None, weights=None):
    """
    Setup the selection calculation
    if select_count [0,1], select that percent of samples
    if select_count >= 1, select that number of samples
    """
    num_vars, num_samples = data["data"].shape
    logging.info("Sample Count %d", num_samples)
    logging.info("Variant Count %d", num_vars)

    select_count = num_samples if select_count < 0 \
                   else max(1, int(num_samples * select_count) if select_count < 1 \
                   else int(select_count))
    logging.info("Selecting %d samples", select_count)

    vcf_samples = data['samples']
    if isinstance(data, h5py.File):
        vcf_samples = data['samples'][:].astype(str)
    else:
        vcf_samples = vcf_samples.astype(str)

    # Build masks
    sample_mask = np.ones(num_samples, dtype='bool')

    if subset:
        sample_mask = np.isin(vcf_samples, subset)
        logging.info("Subsetting to %d of %d samples", sample_mask.sum(), num_samples)
    if exclude:
        logging.info(f"Excluding {len(exclude)} samples")
        sample_mask &= ~np.isin(vcf_samples, exclude)
    if subset or exclude:
        logging.info(f"Ending with {sample_mask.sum()} samples")

    sample_weights = None
    if weights is not None:
        logging.info("Setting %d weights", len(weights))
        sample_weights = np.ones(num_samples)
        for pos, i in enumerate(vcf_samples):
            if i in weights.index:
                sample_weights[pos] = weights.loc[i]

    matrix = data['data']

    if isinstance(data, h5py.File) and is_memsafe(matrix.shape):
        logging.info("Dataset small enough to hold in memory")
        matrix = matrix[:]

    return greedy_select(matrix, data["var_count"][:], select_count, vcf_samples, sample_mask, sample_weights)


def write_append_hdf5(cur_part, out_name, is_first=False, calc_af=False):
    """
    Consolidate parts into the hdf5 file
    """
    n_cols = cur_part['GT'].shape[1]
    # Consider turning 1e6 into a parameter?
    # Or maybe just turning chunking off..
    c_size = (max(1, int(1e6 / 4 / n_cols)), n_cols)

    if is_first:
        with h5py.File(out_name, 'w') as hf:
            hf.create_dataset('samples', data=cur_part["samples"], compression="lzf", chunks=True, maxshape=(None, ))
            if not calc_af:
                hf.create_dataset('data',
                              data=cur_part["GT"],
                              compression="lzf",
                              dtype='bool',
                              chunks=c_size,
                              maxshape=(None, n_cols))
            else:
                hf.create_dataset('data',
                                  data=cur_part["GT"] * cur_part["AF"],
                                  compression="lzf",
                                  dtype='float32',
                                  chunks=c_size,
                                  maxshape=(None, n_cols))
        return

    with h5py.File(out_name, 'a') as hf:
        hf["data"].resize((hf["data"].shape[0] + cur_part["GT"].shape[0]), axis=0)
        if not calc_af:
            hf["data"][-cur_part["GT"].shape[0]:] = cur_part["GT"]
        else:
            hf["data"][-cur_part["GT"].shape[0]:] = cur_part["GT"] * cur_part["AF"]

def add_varcount_to_h5(fn, var_count):
    """
    Add the new var_count dataset to the h5 file
    """
    with h5py.File(fn, 'a') as hf:
        hf.create_dataset('var_count', data=var_count, compression="lzf")

#pylint: disable=too-many-statements
def load_files(in_files, lowmem=None, buffer=32768, calc_af=False):
    """
    Load and concatenate multiple files
    if lowmem is provided, in_files are concatenated into an hdf5. This may be a little slower,
    but will use less memory since it won't need to hold two copies of the
    data in memory during the concatenation step.
    If the af_matrix is going to be used, lowmem will also create that per-in_file and write it to the hdf5 file
    """
    logging.info(f"Loading {len(in_files)} files")
    if lowmem == 1:
        return h5py.File(in_files[0], 'r')

    samples = None
    gt_parts = []
    af_parts = []

    load_row_count = 0
    load_buffer_count = 0

    is_first = True
    var_count = None
    for load_count, i in enumerate(in_files):
        if i.endswith((".vcf.gz", ".vcf")):
            dat = read_vcf(i, lowmem is not None, buffer)
        elif i.endswith(".jl"):
            dat = joblib.load(i)
        else:
            logging.error("Unknown filetype %s. Expected `.vcf[.gz]`, `.jl`", i)
            sys.exit(1)

        if samples is None:
            samples = dat['samples'].astype('S')

        upack = np.unpackbits(dat['GT'], axis=1, count=len(dat['samples'])).astype(bool)
        uninf_filter = upack.any(axis=1)
        logging.debug("fitering %d uninformative variants", (~uninf_filter).sum())

        gt_parts.append(upack[uninf_filter])
        af_parts.append(dat['AF'][uninf_filter])
        if var_count is None:
            var_count = gt_parts[-1].sum(axis=0)
        else:
            var_count += gt_parts[-1].sum(axis=0)

        # lowmem buffering
        m_count = dat["GT"].shape[0]
        load_row_count += m_count
        load_buffer_count += m_count
        if lowmem is not None and (load_buffer_count >= buffer or load_count + 1 == len(in_files)):
            logging.debug("dumping chunk with %d vars", load_buffer_count)
            cur_chunk = None
            if len(gt_parts) > 1:
                cur_chunk = {'GT': np.concatenate(gt_parts), 'samples': samples, 'AF': np.concatenate(af_parts)}
            else:
                cur_chunk = {'GT': gt_parts[0], 'samples': samples, 'AF': af_parts[0]}

            write_append_hdf5(cur_chunk, lowmem, is_first, calc_af)

            load_buffer_count = 0
            is_first = False
            gt_parts = []
            af_parts = []

        logging.debug("Loaded %d of %d (%.2f%%) with %d vars (%d buff)", load_count + 1, len(in_files),
                      (load_count + 1) / len(in_files) * 100, load_row_count, load_buffer_count)
        logging.debug("Average of %d variants per-sample", np.mean(var_count))

    if lowmem is not None:
        # Need to write the var_count to the h5file
        add_varcount_to_h5(lowmem, var_count)
        return h5py.File(lowmem, 'r')

    ret = {'samples': samples}
    ret['data'] = np.concatenate(gt_parts) if len(gt_parts) > 1 else gt_parts[0]
    ret['var_count'] = var_count
    if calc_af:
        logging.info("Calculating AF Matrix")
        af_arr = np.concatenate(af_parts) if len(af_parts) > 1 else af_parts[0]
        ret["data"] = ret["data"] * af_arr

    return ret
#pylint: enable=too-many-statements

###################
# Input utilities #
###################
def parse_sample_lists(argument):
    """
    Parse the --exclude/--subset arguments
    """
    ret = []
    if not argument:
        return ret
    for i in argument:
        if os.path.exists(i):
            with open(i, 'r') as fh:
                ret.extend([_.strip() for _ in fh])
        else:
            ret.extend(i.split(","))
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
    parser = argparse.ArgumentParser(prog="select",
                                     description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_files", nargs="*", type=str, help="Input VCF or jl files")
    parser.add_argument("-c",
                        "--count",
                        type=float,
                        default=0.02,
                        help="Number of samples to select as a percent if <1 or count if >=1 or -1 for all (%(default)s)")
    parser.add_argument("-o", "--out", type=str, default="/dev/stdout", help="Output file (stdout)")
    parser.add_argument("--debug", action="store_true", help="Verbose logging")

    scoreg = parser.add_argument_group("Scoring Arguments")
    scoreg.add_argument("--af", action="store_true", help="Weigh variants by allele frequency")
    scoreg.add_argument("--weights", type=str, default=None, help="Tab-delimited file of sample weights")
    scoreg.add_argument("--subset",
                        type=str,
                        default=None,
                        action='append',
                        help="Filename with or Comma-separated list of samples to analyze")
    scoreg.add_argument("--exclude",
                        type=str,
                        default=None,
                        action='append',
                        help="Filename with or Comma-separated list of samples to exclude selection")

    mperfg = parser.add_argument_group("Memory Arguments")
    mperfg.add_argument("--lowmem",
                        type=str,
                        default=None,
                        help="Name of concatenated hdf5 file to create/use (%(default)s)")
    mperfg.add_argument("--buffer",
                        type=int,
                        default=32768,
                        help="Number of variants to buffer during concatenation (%(default)s)")
    mperfg.add_argument("--maxmem",
                        type=int,
                        default=2,
                        help="Maximum amount of memory in (GB). 0 keeps data in hdf5 (%(default)s)")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    # Validate inputs
    if [_ for _ in args.in_files if _.endswith(".hdf5")] and len(args.in_files) > 1:
        logging.error("Cannot provide hdf5 with multiple input files")
        sys.exit(1)

    if len(args.in_files) == 0:
        if not args.lowmem:
            logging.error("No input files provided")
            sys.exit(1)
        args.in_files = [args.lowmem]
        args.lowmem = 1

    if len(args.in_files) == 1 and args.in_files[0].endswith(".hdf5") and not args.lowmem:
        logging.info("Switching on lowmem for hdf5 input")
        args.lowmem = 1

    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))
    return args


def select_main(cmdargs):
    """
    Main
    """
    global MAXMEM  # pylint: disable=global-statement
    args = parse_args(cmdargs)

    data = load_files(args.in_files, args.lowmem, args.buffer, args.af)
    if data['data'].dtype == bool and args.af:
        logging.critical("HDF5 file doesn't appear to be created with --af weighted scores, remove --af or recreate hdf5")
        sys.exit(1)
    if data['data'].dtype != bool and not args.af:
        logging.critical("HDF5 file appears to be created with --af weighted scores, add --af or recreate hdf5")

    args.subset = parse_sample_lists(args.subset)
    args.exclude = parse_sample_lists(args.exclude)
    args.weights = parse_weights(args.weights)

    MAXMEM = args.maxmem
    with open(args.out, 'w') as fout:
        fout.write("sample\tvar_count\tnew_count\ttot_captured\tpct_captured\n")
        m_iter = run_selection(data, args.count, args.subset, args.exclude, args.weights)
        for result in m_iter:
            logging.info("Selected %s (%.1f%% of variants)", result[0], result[4] * 100)
            fout.write("\t".join([str(_) for _ in result]) + '\n')
            fout.flush()

    logging.info("Finished utmos")
