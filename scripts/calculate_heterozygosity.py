"""
calculate_heterozygosity.py
Daniel Cotter

calculates heterozygosity at every locus in a vcf file
requires a vcf that is subset for samples of interest
"""
# import required modules
import sys
import csv
from os import path
import argparse
import vcf  # installed via conda as pyvcf


# define useful functions
def heterozygosity_site(counts, total):

    return freq


# parse command line arguments
parser = argparse.ArgumentParser(description="Calculates heterozygosity " +
                                 "at every locus in a vcf file")

parser.add_argument("--vcf", nargs='?', default=sys.stdin,
                    help="Input VCF file. Can be gzipped. Reads from stdin" +
                    " by default.")
parser.add_argument("--chrom", default=None,
                    help="Provide the name of the scaffold or chromosome to" +
                    " be analyzed. Otherwise, will output all sites.")
parser.add_argument("--outuput", nargs='?', default=sys.stdout,
                    help="Output file location. Default is stdout.")

# Print help/usage if no arguments are supplied
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

# begin script ----------------------------------------------------------------

# initialize the vcf reader object
vcf_reader = vcf.Reader(filename=args.vcf)
