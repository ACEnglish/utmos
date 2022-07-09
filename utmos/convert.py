import sys
import logging
import argparse

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
    parser.add_argument("-c", "--compress", type=int, default=5,
                        help="joblib compress level 1-9 (%(default)s)")
    # parser.add_argument --af pull or calculate allele frequencies
    args = parser.parse_args(args)
    truvari.setup_logging()
    return args

def read_vcf(in_file):
    """
    Read a vcf's genotypes and return numpy arrays
    """
    logging.info("Reading VCF")
    data = allel.read_vcf(in_file, fields=["calldata/GT", "samples"])
    logging.info(f"Converting genotypes to bool {data['calldata/GT'].shape}")
    gts = allel.GenotypeArray(data["calldata/GT"])
    is_het = gts.is_het()
    logging.info(f"{is_het.sum().sum()} hets")
    is_hom = gts.is_hom_alt()
    logging.info(f"{is_hom.sum().sum()} homs")
    v_count = is_het | is_hom
    del(data["calldata/GT"])
    data["GT"] = v_count
    return data

def cvt_main(cmdargs):
    """
    Main 
    """
    args = parse_args(cmdargs)
    truvari.setup_logging()
    # fields = optionaly AF (make sure its consistent with Number=[./1 whatever]
    data = read_vcf(args.in_file)
    logging.info("Saving genotypes")
    joblib.dump(data, args.out_file, compress=args.compress)
    logging.info("Finished conversion")
