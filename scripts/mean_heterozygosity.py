"""
calculate_heterozygosity.py
Daniel Cotter

calculates heterozygosity at every locus in a vcf file
requires a vcf that is subset for samples of interest
"""
# import required modules
import sys
import csv
import argparse
import re
import os
from scipy import stats
import numpy as np


# define useful functions -----------------------------------------------------
def get_pop_chr(filepath):
    """
    Return the population code and chromosome code from the filename provided
    """
    pop_pattern = re.compile(r'[A-Z]{3}')
    chr_pattern = re.compile(r'chr([XY]{1}|\d+)')
    filename = os.path.basename(filepath)
    return [pop_pattern.search(filename)[0],
            chr_pattern.search(filename)[0]]


# parse command line arguments
parser = argparse.ArgumentParser(description="Determines mean heterozygosity" +
                                             " for all input files")
parser.add_argument("--input_files", nargs='+',
                    help="List of all input files for means to be calculated.")
parser.add_argument("--output", nargs='?', required=True,
                    help="Merged output file.")

# Print help/usage if no arguments are supplied
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

# begin script ----------------------------------------------------------------
# loop through all provided files and write each line to output
with open(args.output, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')

    for file in args.input_files:
        pop, chrom = get_pop_chr(file)
        # open the file and grab all heterozygosity values
        with open(file, 'r') as f:
            het_vals = []
            for line in f:
                line = line.strip().split('\t')
                het_vals.append(float(line[3]))
            n = len(het_vals)
            avg = np.mean(het_vals)
            sem = stats.sem(het_vals)
            h = sem * stats.t.ppf((1 + 0.95) / 2, n - 1)
            writer.writerow([pop, chrom, avg,
                             avg - h, avg + h, avg - sem, avg + sem])

        print("{} merged".format(file))
