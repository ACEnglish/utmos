"""
convert vcfs to lightweight numpy arrays
"""
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
    parser.add_argument("-p", "--nopackbits", action="store_false",
                        help="Use numpy.packbits to make output smaller (%(default)s)")
    parser.add_argument("--af", action="store_true",
                        help="Calcluate allele frequencies (%(default)s)")

    args = parser.parse_args(args)
    truvari.setup_logging()
    return args

def read_vcf(in_file, lowmem=False, allele_freq=False, packbits=False):
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
    num_hets = is_het.sum().sum()
    logging.info(f"{num_hets} hets")

    is_hom = gts.is_hom_alt()
    num_homs = is_hom.sum().sum()
    logging.info(f"{num_homs} homs")
    v_count = is_het | is_hom

    if lowmem:
        data = {"GT": v_count, "samples": data["samples"][:].astype(str)}
    else:
        del data["calldata/GT"]
        data["GT"] = v_count

    if packbits:
        data["GT"] = np.packbits(data["GT"], axis=1)
        data['packedbits'] = True
    else:
        data['packedbits'] = False

    af = None
    if allele_freq:
        af = gts.count_alleles().to_frequencies()[:, 1]

    if allele_freq:
        data["AF"] = af

    data["stats"] = {'num_het': num_hets, 'num_hom': num_homs}
    return data

def cvt_main(cmdargs):
    """
    Main
    """
    args = parse_args(cmdargs)
    # fields = optionaly AF (make sure its consistent with Number=[./1 whatever]
    data = read_vcf(args.in_file, args.lowmem, args.af, args.nopackbits)
    logging.info("Saving genotypes")
    joblib.dump(data, args.out_file, compress=args.compress)
    logging.info("Finished conversion")
