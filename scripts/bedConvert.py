"""
bedConvert.py
Daniel Cotter

This script converts single position diversity coordinates into
BED format for easy filtering with bedtools.
Last Updated: 5/7/2018
"""

# Import required modules
from sys import argv
import csv

# parse argument variables
script, input_file, output_file = argv

# read in diversity files
with open(input_file, 'rU') as f:
    positions = list(csv.reader(f, delimiter='\t'))

# modify coordinates and write out diversity files
with open(output_file, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t', newline='\n')
    for info in positions:
        writer.writerow(['chr' + info[0], str(int(info[1]) - 1),
                         info[1], info[2]])
