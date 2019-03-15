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


# define useful functions -----------------------------------------------------
def get_pop_chr(filepath):
    """
    Return the population code and chromosome code from the filename provided
    """
    pop_pattern = re.compile(r'[A-Z]{3}')
    chr_pattern = re.compile(r'chr([XY]{1}|\d+)')
    filename = os.path.basename(filepath)
    return [pop_pattern.match(filename), chr_pattern.match(filename)]


# parse command line arguments
parser = argparse.ArgumentParser(description="Determines mean heterozygosity" +
                                             " for all input files")
parser.add_argument("--input_files", nargs='+',
                    help="List of all input files for means to be calculated.")
parser.add_argument("--output", nargs='?', default=True,
                    help="Merged output file. Default is stdout.")

# Print help/usage if no arguments are supplied
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

# begin script ----------------------------------------------------------------
# initialize data dictionary
data = {}

# loop through all provided files
for file in args.input_files:
    pop, chrom = get_pop_chr(file)
    if pop not in data:
        data[pop] = {}
    # open the file and grab all heterozygosity values
    with open(file, 'r') as f:
        het_vals = []
        for line in f:
            line = line.strip().split('\t')
            het_vals.append(float(line[3]))
    # write the pop, chr, and het values to a list in the dictionary
    data[pop][chrom] = [pop, chr, ','.join(het_vals)]

# reformat data to be output
results = []
for pop in data:
    for chr in pop:
        results.append(chr)

# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t')
    for row in results:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        for row in results:
            writer.writerow(row)
