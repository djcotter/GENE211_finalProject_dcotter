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
SEX = ['individuals']
POPS = SUBPOPULATIONS
CHROM = ['chr8']

# Rules -----------------------------------------------------------------------

# Rule ALL
rule all:
    input:
        expand(
            path.join('results',
                      '{chr}_merged_heterozygosity_{sex}_{filter_iter}.txt'),
            chr=['chr8', 'chrX', 'chrY'], sex=SEX, filter_iter=FILTER)

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
# an example of python3 urllib
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
        path.join('data', 'subset_{chr}_{pop}_{sex}.vcf')
    params:
        script = path.join('scripts', 'calculate_heterozygosity.py'),
        chrom = lambda wildcards: wildcards.chr[3:]
    output:
        temp(path.join('results',
                       '{chr}_{pop}_{sex}_heterozygosity_by_site.txt'))
    conda:
        path.join('envs', 'calculate_diversity.yml')
    shell:
        "python {params.script} --vcf {input} --chrom {params.chrom} "
        "--output {output}"

# convert output of diversity script to BED format
rule convert_to_BED:
    input:
        path.join('results', '{chr}_{pop}_{sex}_heterozygosity_by_site.txt')
    params:
        script = path.join('scripts', 'bedConvert.py')
    output:
        path.join('results', '{chr}_{pop}_{sex}_heterozygosity_by_site.bed')
    shell:
        "python {params.script} {input} {output}"

# use UCSC interval files to create desired filter
rule create_filter:
    input:
        lambda wildcards: expand(
            path.join('data', 'raw_filters', '{filter_type}'),
            filter_type=config["filters"][wildcards.filter_iter])
    params:
        chrom = lambda wildcards: wildcards.chr
    output:
        path.join('filters', '{chr}_complete_{filter_iter}.bed')
    conda:
        path.join('envs', 'calculate_diversity.yml')
    shell:
        "cat {input} | grep {params.chrom} | sort -k1,1 -k2,2n | "
        "cut -f1,2,3 | bedtools merge -i stdin > {output}"

# filter the BED output of the diversity file
rule filter_heterozygosity_by_site:
    input:
        het_data = path.join('results', '{chr}_{pop}_{sex}' +
                             '_heterozygosity_by_site.bed'),
        filter = path.join('filters', '{chr}_complete_{filter_iter}.bed')
    output:
        temp(path.join('results', '{chr}_{pop}_{sex}' +
                       '_heterozygosity_by_site_{filter_iter}.bed'))
    conda:
        path.join('envs', 'calculate_diversity.yml')
    shell:
        "bedtools subtract -a {input.het_data} -b {input.filter} > {output}"

# given a file of heterozygosity_by_site return a single mean with 95% CI
rule merge_heterozygosity_output:
    input:
        lambda wildcards: expand(
            path.join('results', '{chr}_{pop}_{sex}' +
                      '_heterozygosity_by_site_{filter_iter}.bed'),
            chr=wildcard.chr, pop=POPS, sex=wildcards.sex,
            filter_iter=wildcards.filter_iter)
    params:
        script = path.join('scripts', 'merge_heterozygosity.py'),
    output:
        path.join('results',
                  '{chr}_merged_heterozygosity_{sex}_{filter_iter}.txt')
    shell:
        "python {params.script} --input_files {input.autosomes} "
        "{input.chrX} {input.chrY} --output {output}"
