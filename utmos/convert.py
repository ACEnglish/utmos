"""
convert vcfs to lightweight numpy arrays
"""
import json
import logging
import argparse
import tempfile

import h5py
import allel
import joblib
import truvari
import numpy as np


def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="convert",
                                     description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_file", type=str, help="Input VCF")
    parser.add_argument("out_file", type=str, help="Output joblib")
    parser.add_argument("--no-singleton", action="store_true",
                        help="Remove singleton variants")
    parser.add_argument("--lowmem",
                        action="store_true",
                        help="Lower memory usage with hdf5 temporary files (%(default)s)")
    parser.add_argument("-B",
                        "--buffer",
                        type=int,
                        default=50000,
                        help="Number of variants read at a time by scikit-allel (%(default)s)")
    parser.add_argument("-c", "--compress", type=int, default=5, help="joblib compress level 1-9 (%(default)s)")

    args = parser.parse_args(args)
    truvari.setup_logging()
    logging.info("Params:\n%s", json.dumps(vars(args), indent=4))
    return args


def read_vcf(in_file, lowmem=False, chunk_length=2000, no_singleton=False):
    """
    Read a vcf's genotypes and return numpy arrays
    """
    logging.info("Reading VCF")
    if lowmem:
        tmpfile = tempfile.mkstemp()[1]
        allel.vcf_to_hdf5(in_file, tmpfile, fields=["calldata/GT", "samples"], chunk_length=chunk_length)
        data = h5py.File(tmpfile)
    else:
        data = allel.read_vcf(in_file, fields=["calldata/GT", "samples"], chunk_length=chunk_length)

    logging.info(f"Converting genotypes to bool {data['calldata/GT'].shape}")

    gts = allel.GenotypeArray(data["calldata/GT"])
    if no_singleton:
        ac = gts.count_alleles()
        sing = ac.is_singleton(allele=1) | ac.is_singleton(allele=0)
        logging.info("Removing %d singletons", sing.sum())
        gts = gts[~sing]

    is_het = gts.is_het()
    num_hets = is_het.sum().sum()
    logging.info(f"{num_hets} hets")

    is_hom = gts.is_hom_alt()
    num_homs = is_hom.sum().sum()
    logging.info(f"{num_homs} homs")
    v_count = is_het | is_hom

    logging.info("Calculating AFs")
    # Use maximum non-reference allele frequency
    af = gts.count_alleles().to_frequencies()[:, 1:].max(axis=1)
    # Needs to be reshaped for future multiplications
    af = af.reshape(af.shape[0], 1)

    if lowmem:
        data = {"GT": v_count, "samples": data["samples"][:].astype(str)}
    else:
        del data["calldata/GT"]
        data["GT"] = v_count
    data["AF"] = af
    data["GT"] = np.packbits(data["GT"], axis=1)

    data["stats"] = {'num_het': num_hets, 'num_hom': num_homs}
    return data


def cvt_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    data = read_vcf(args.in_file, args.lowmem, args.buffer, args.no_singleton)
    logging.info("Saving genotypes")
    joblib.dump(data, args.out_file, compress=args.compress)
    logging.info("Finished conversion")
