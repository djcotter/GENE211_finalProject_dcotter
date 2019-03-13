"""
population_parser.py
Daniel Cotter

Take a reference panel from 1000 genomes and parse by pops and gender
---------------------------------------------------------------------

usage is population_parser.py input_file output_directory
"""

# import statements

from sys import argv
import csv
from os import path

# -----------------------------------------------------------------------------
# declare argument variables
script, input_file, output_directory = argv

# open the input population panel file
with open(input_file, 'rU') as f:
    samples = list(csv.reader(f, delimiter='\t'))

# initialize lists for male and female individuals
males = []
females = []

# loop through all samples in the panel file and assign them male or female
for samp in samples:
    if samp[3] == 'male':
        males.append(samp[0])
    elif samp[3] == 'female':
        females.append(samp[0])

# create dictionaries with seperate lists of individuals by pop
pop_dict_males = {}
pop_dict_females = {}

# loop through samples and identify population codes
for samp in samples:
    if samp[3] == 'male':
        if samp[2] in pop_dict_males:
            pop_dict_males[samp[2]].append(samp[0])
        else:
            pop_dict_males[samp[2]] = [samp[0]]
    elif samp[3] == 'female':
        if samp[2] in pop_dict_females:
            pop_dict_females[samp[2]].append(samp[0])
        else:
            pop_dict_females[samp[2]] = [samp[0]]

# create dictionaries with seperate lists of individuals by subpop
subpop_dict_males = {}
subpop_dict_females = {}

# repeat the above using the population code and not the superpopulation code
for samp in samples:
    if samp[3] == 'male':
        if samp[1] in subpop_dict_males:
            subpop_dict_males[samp[1]].append(samp[0])
        else:
            subpop_dict_males[samp[1]] = [samp[0]]
    elif samp[3] == 'female':
        if samp[1] in subpop_dict_females:
            subpop_dict_females[samp[1]].append(samp[0])
        else:
            subpop_dict_females[samp[1]] = [samp[0]]

# output the files as lists of individuals
# loop through files and print parsed lists of individuals
for pop in pop_dict_females:
    temp_path = path.join(output_directory, pop + '_females')
    with open(temp_path, 'w') as file:
        for samp in pop_dict_females[pop]:
            file.write(samp + '\n')
    temp_path2 = path.join(output_directory, pop + '_males')
    with open(temp_path2, 'w') as file:
        for samp in pop_dict_males[pop]:
            file.write(samp + '\n')
    temp_path3 = path.join(output_directory, pop + '_individuals')
    with open(temp_path3, 'w') as file:
        for samp in pop_dict_females[pop]:
            file.write(samp + '\n')
        for samp in pop_dict_males[pop]:
            file.write(samp + '\n')

# loop through the files and print parsed lists of individuals
for pop in subpop_dict_females:
    temp_path = path.join(output_directory, pop + '_females')
    with open(temp_path, 'w') as file:
        for samp in subpop_dict_females[pop]:
            file.write(samp + '\n')
    temp_path2 = path.join(output_directory, pop + '_males')
    with open(temp_path2, 'w') as file:
        for samp in subpop_dict_males[pop]:
            file.write(samp + '\n')
    temp_path3 = path.join(output_directory, pop + '_individuals')
    with open(temp_path3, 'w') as file:
        for samp in subpop_dict_females[pop]:
            file.write(samp + '\n')
        for samp in subpop_dict_males[pop]:
            file.write(samp + '\n')

# write a file with a list of ALL individuals split by sex
temp_path = path.join(output_directory, 'ALL' + '_females')
with open(temp_path, 'w') as file:
    for samp in females:
        file.write(samp + '\n')

temp_path = path.join(output_directory, 'ALL' + '_males')
with open(temp_path, 'w') as file:
    for samp in males:
        file.write(samp + '\n')

temp_path = path.join(output_directory, 'ALL' + '_individuals')
with open(temp_path, 'w') as file:
    for samp in females:
        file.write(samp + '\n')
    for samp in males:
        file.write(samp + '\n')
