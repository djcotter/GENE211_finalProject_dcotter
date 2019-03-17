"""
calculate_heterozygosity.py
Daniel Cotter

calculates heterozygosity at every locus in a vcf file
requires a vcf that is subset for samples of interest
"""
# import required modules -----------------------------------------------------
import sys
import csv
import argparse
import vcf  # installed via conda as pyvcf
import numpy as np
import re


# define useful functions -----------------------------------------------------
def expected_heterozygosity(alt_counts, total_alleles):
    """
    Given a list of the counts of each alternate allele and the total number
    return the expected heterozygosity at a site.
    """
    alt_sum = np.sum(alt_counts)
    all_counts = [total_alleles - alt_sum] + alt_counts
    freqs = [float(x / total_alleles) for x in all_counts]
    h = 1 - np.sum(np.square(freqs))
    return h


# parse command line arguments
parser = argparse.ArgumentParser(description="Calculates heterozygosity " +
                                 "at every locus in a vcf file")

parser.add_argument("--vcf", nargs='?', default=sys.stdin,
                    help="Input VCF file. Can be gzipped. Reads from stdin" +
                    " by default.")
parser.add_argument("--chrom", default=None,
                    help="Provide the name of the scaffold or chromosome to" +
                    " be analyzed. Otherwise, will output all sites.")
parser.add_argument("--output", nargs='?', default=True,
                    help="Output file location. Default is stdout.")

# Print help/usage if no arguments are supplied
if len(sys.argv) == 1:
    parser.print_usage()
    sys.exit(1)

args = parser.parse_args()

# begin script ----------------------------------------------------------------

# initialize the vcf reader object
vcf_reader = vcf.Reader(filename=args.vcf)

# initialize results array
results = []
chr_num = re.compile(r'(^\d{1,2}|^[XY])')
for record in vcf_reader:
    # test that the chromosome is correct
    if record.CHROM == args.chrom or args.chrom is None:
        # use the info stored in the VCF to calculate heterozygosity
        het = expected_heterozygosity(record.INFO['AC'], record.INFO['AN'])
        chrom_string = 'chr' + chr_num.match(record.CHROM).group(1)
        results.append([chrom_string, record.POS, het])

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
