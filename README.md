# **GENE211 Final Project**
**Daniel Cotter**
**17 March 2019**

This code is associated with a final project for GENE211 at Stanford University. The files here can be run using `Snakemake`

## Instructions
1. Make sure you have `conda` installed because the script depends on specific environments. You can install `conda` [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

2. Install `snakemake`. The easiest way to do this is using `conda`:
```bash
conda install snakemake
```
3. Run the pipeline by changing into the project directory and running:
```bash
snakemake --use-conda
```

### Notes
`conda` has been running into some issues with `R` dependencies that might make the last few jobs crash.
