"""
1000 Genomes Heterozygosity
Daniel Cotter

analyze average heterozygosity across 22 autosomes
---------------------------------------------------------------------

Requires:
    conda
        to activate the included mut_recom.yml environment
    python3
        snakemake needs python >= 3.5 to run
"""

# Import required packages ----------------------------------------------------

import json
from os import path


# Global Configurations -------------------------------------------------------

# path to configuration file
configfile: 'config.yml'

# parse the populations from populations.json
POPULATIONS = sorted(json.load(open(config['POP_CODES']))['Populations'])
SUBPOPULATIONS = sorted(json.load(open(config['POP_CODES']))['Subpopulations'])

# specify filter, window size, sex, and populations to be analyzed
FILTER = ['filter1']
WINDOW = ['100kb']
SEX = ['individuals']
POPS = ['YRI']
CORRECTION = ['uncorrected']
CHROM = ['22']


# Rules -----------------------------------------------------------------------

# Rule ALL
rule all:
    input:
        expand(path.join('04_window_analysis', 'results',
                         '{pop}_{sex}_{chr}_{filter_iter}_{window}' +
                         '_diversity.bed'),
               pop=POPS, sex=SEX, chr=CHROM,
               filter_iter=FILTER, window=WINDOW)

# download VCF files
rule download_VCF_files:
    params:
        url = lambda wildcards:
            config['chromosomes'][wildcards.chr]['data_url']
    output:
        path.join('data', '1000genomes_{chr}.vcf.gz')
    run:
        import urllib
        urllib.request.urlretrieve(params[0], output[0])

# download the panel file
rule download_panel_file:
    params:
        url = config['panel']['data_url']
    output:
        config['panel']['project_path']
    run:
        import urllib
        urllib.request.urlretrieve(params[0], output[0])

# parse the populations from the panel file
rule parse_populations:
    input:
        config['panel']['project_path']
    params:
        script = path.join('scripts', 'population_parser.py'),
        out_dir = path.join('populations')
    output:
        out_pops = expand(path.join('populations', '{pops}_{group}'),
                          pops=POPULATIONS + SUBPOPULATIONS,
                          group=['males', 'females', 'individuals'])
    shell:
        "python {params.script} {input} {params.out_dir}"

# split the input VCFs by population
rule subset_and_filter_VCF:
    input:
        population = path.join('populations', '{pop}_{sex}'),
        vcf = path.join('data', '1000genomes_{chr}.vcf.gz')
    output:
        temp(path.join('data', 'subset_{chr}_{pop}_{sex}.vcf'))
    conda:
        path.join('envs', 'calculate_diversity.yml')
    shell:
        "bcftools view -S {input.population} {input.vcf} | "
        "bcftools view -Ou -m2 -M2 -v snps | bcftools view "
        "-Ov --min-ac 1:minor > {output}"

# calculate heterozygostiy on the autosomes by site
rule heterozygosity_by_site:
    input:
        population = path.join('populations', '{pop}_{sex}'),
        vcf = path.join('data', 'subset_{chr}_{pop}_{sex}.vcf')
    params:
        calc_pi = path.join('scripts', 'calculate_heterozygosity.py'),
        out_dir = path.join('02_diversity_by_site', 'results/'),
        chrom = lambda wildcards: wildcards.chr[3:]
    output:
        temp(path.join('results',
                       '{pop}_{sex}_{chr}_heterozygosity_by_site.txt'))
    conda:
        path.join('envs', 'calculate_diversity.yml')
    shell:
        "python {params.calc_pi} --vcf {input.vcf} "
        "--population_lists {input.population} --chrom_inc {params.chrom} "
        "--out_directory {params.out_dir}"

# convert output of diversity script to BED format
rule convert_to_BED:
    input:
        path.join('02_diversity_by_site', 'results',
                  '{pop}_{sex}_{chr}_pi_output_by_site.txt')
    params:
        script = path.join('scripts', 'bedConvert.py')
    output:
        path.join('02_diversity_by_site', 'results',
                  '{pop}_{sex}_{chr}_pi_output_by_site.bed')
    shell:
        "python {params.script} {input} {output}"

# use UCSC interval files to create desired filter
rule create_filter:
    input:
        lambda wildcards: expand(
            '03_filters/raw_filters/{filter_type}',
            filter_type=config["filters"][wildcards.filter_iter])
    params:
        chrom = lambda wildcards: wildcards.chr
    output:
        path.join('03_filters', 'results',
                  'complete_{chr}_{filter_iter}.bed')
    shell:
        "cat {input} | grep {params.chrom} | sort -k1,1 -k2,2n | "
        "cut -f1,2,3 | bedtools merge -i stdin > {output}"

# filter the BED output of the diversity file
rule filter_diversity_by_site:
    input:
        diversity_by_site = path.join('02_diversity_by_site', 'results',
                                      '{pop}_{sex}_{chr}' +
                                      '_pi_output_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_{filter_iter}.bed')
    output:
        path.join('04_window_analysis', 'inputs',
                  '{pop}_{sex}_{chr}_{filter_iter}' +
                  '_pi_by_site.bed')
    shell:
        "bedtools intersect -a {input.diversity_by_site} "
        "-b {input.filtered_callable} > {output}"

rule window_analysis:
    input:
        filtered_diversity = path.join('04_window_analysis', 'inputs',
                                       '{pop}_{sex}_{chr}_{filter_iter}' +
                                       '_pi_by_site.bed'),
        filtered_callable = path.join('04_window_analysis', 'inputs',
                                      'callable_sites_{chr}_' +
                                      '{filter_iter}.bed'),
        windows = path.join('04_window_analysis', 'inputs',
                            '{chr}_{window}_window.bed')
    params:
        window_calcs = path.join('04_window_analysis', 'scripts',
                                 'window_calculations.py'),
        slide = lambda wildcards: '' if \
            config["windows"][wildcards.window]["overlap"] is False else \
            '--sliding ',
        replicates = 1000
    output:
        temp(path.join('04_window_analysis', 'results',
                       '{pop}_{sex}_{chr}_{filter_iter}_{window}' +
                       '_diversity.unfiltered.bed'))
    shell:
        "python {params.window_calcs} --diversity {input.filtered_diversity} "
        "--callable {input.filtered_callable} --windows {input.windows} "
        "{params.slide}--replicates {params.replicates} --output {output}"

rule filter_windows_by_callable_sites:
    input:
        path.join('04_window_analysis', 'results',
                  '{pop}_{sex}_{chr}_{filter_iter}_{window}' +
                  '_diversity.unfiltered.bed')
    params:
        script = path.join('04_window_analysis', 'scripts',
                           'filter_windows_byCallableSites.py'),
        winSize = lambda wildcards:
            config["windows"][wildcards.window]["win_size"]
    output:
        path.join('04_window_analysis', 'results',
                  '{pop}_{sex}_{chr}_{filter_iter}_{window}' +
                  '_diversity.bed')
    shell:
        "python {params.script} --input {input} --windowSize {params.winSize} "
        "--filter 0.1 --output {output}"
