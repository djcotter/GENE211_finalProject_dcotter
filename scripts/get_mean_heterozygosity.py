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
import numpy as np
import re
import os


# define useful functions -----------------------------------------------------
def get_pop_chr(filepath):
    """
    Return the population code and chromosome code from the filename provided
    """
    pop_pattern = re.compile(r'[A-Z]{3}')
    chr_pattern = re.compile(r'chrX|chrY|autosomes')
    filename = os.path.basename(filepath)
    return [pop_pattern.match(filename), chr_pattern.match(filename)]


# Bootstrap functions
def rand_samples(data_array):
    """
    Create a randomly dsitributed set of data based on a given array and
    the length of the given data array
    """
    array_len = len(data_array)
    if array_len > 0:
        indices = np.random.randint(array_len - 1, size=array_len)
        return [data_array[i] for i in indices]
    else:
        return []


def bootstrap_CI_mean(data_array, replicates):
    """
    Bootstraps the data to get a 95% confidence interval of the mean
    returns a list with the mean
    """
    resamples = []
    for i in range(replicates):
        samples = rand_samples(data_array)
        if samples:
            resamples.append(np.mean(samples))
    return [np.nanpercentile(resamples, 2.5),
            np.nanpercentile(resamples, 97.5)]


# parse command line arguments
parser = argparse.ArgumentParser(description="Determines mean heterozygosity" +
                                             " for all input files")
parser.add_argument("--input_files", nargs='+',
                    help="List of all input files for means to be calculated.")
parser.add_argument("--output", nargs='?', default=True,
                    help="Merged output file. Default is stdout.")
parser.add_argument("--replicates", nargs='?', type=int, default=1000,
                    help="number of times the bootstrap " +
                    "should resample. Default is 1000.")

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
        for line in file:
            het_vals.append(float(line[2]))
    # write the mean and CI to a line in the dictionary
    data[pop][chrom] = [pop, chr, np.mean(het_vals)] + \
        bootstrap_CI_mean[het_vals, args.replicates]

# reformat data to be output
results = []
for pop in data:
    for chr in pop:
        results.append(chr)

# write the results to output_file or standard out depending on args
if args.output is True:
    writer = csv.writer(sys.stdout, delimiter='\t', newline='\n')
    for row in results:
        writer.writerow(row)
else:
    with open(args.output, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', newline='\n')
        for row in results:
            writer.writerow(row)
