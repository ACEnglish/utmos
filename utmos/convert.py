import sys
import logging
import argparse
import tempfile

import h5py
import allel
import joblib
import truvari

def parse_args(args):
    """
    Pull the command line parameters
    """
    parser = argparse.ArgumentParser(prog="convert", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_file", type=str,
                        help="Input VCF")
    parser.add_argument("out_file", type=str,
                        help="Output joblib")
    parser.add_argument("--lowmem", action="store_true",
                        help="Lower memory usage with hdf5 temporary files (%(default)s)")
    parser.add_argument("-c", "--compress", type=int, default=5,
                        help="joblib compress level 1-9 (%(default)s)")
    # parser.add_argument --af pull or calculate allele frequencies
    args = parser.parse_args(args)
    truvari.setup_logging()
    return args

def read_vcf(in_file, lowmem=False):
    """
    Read a vcf's genotypes and return numpy arrays
    """
    logging.info("Reading VCF")
    if lowmem:
        tmpfile = tempfile.mkstemp()[1]
        allel.vcf_to_hdf5(in_file, tmpfile, fields=["calldata/GT", "samples"])
        data = h5py.File(tmpfile)
    else:
        data = allel.read_vcf(in_file, fields=["calldata/GT", "samples"])

    logging.info(f"Converting genotypes to bool {data['calldata/GT'].shape}")
    gts = allel.GenotypeArray(data["calldata/GT"])
    is_het = gts.is_het()
    logging.info(f"{is_het.sum().sum()} hets")
    is_hom = gts.is_hom_alt()
    logging.info(f"{is_hom.sum().sum()} homs")
    v_count = is_het | is_hom

    if lowmem:
        data = {"GT": v_count, "samples": data["samples"][:]}
    else:
        del(data["calldata/GT"])
    data["GT"] = v_count
    return data

def cvt_main(cmdargs):
    """
    Main 
    """
    args = parse_args(cmdargs)
    # fields = optionaly AF (make sure its consistent with Number=[./1 whatever]
    data = read_vcf(args.in_file, args.lowmem)
    logging.info("Saving genotypes")
    joblib.dump(data, args.out_file, compress=args.compress)
    logging.info("Finished conversion")
