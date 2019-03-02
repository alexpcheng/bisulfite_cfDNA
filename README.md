# bisulfite_cfDNA
Custom scripts used for data analysis of microbial and human derived cell-free DNA after bisulfite treatment published in [...]
Major data analysis pipeline is written using the snakemake workflow (https://snakemake.readthedocs.io/en/stable/).

## Required pipelines
These pipelines are called in custom scripts. They should be obtained and tested before running entire pipeline.
- Bismark
- CrossMap
- Methpipe
- bwa-meth
- MethylDackel

## Folder structure
Most scripts are run relative to the current directory, so initialize your workspace similarly to the repo

## Main pipeline
'''
cd [workspace]
source activate [conda environment]
snakemake
'''
## Methylation references pipeline
'''
cd Methylation_References_snaked
source activate [conda environment]
snakemake
'''
## Run scripts for specific tasks
'''
cd [workspace]/Bin

# Donor fractions
bash donor_fractions_par.sh [make sure you have enough CPUs available, or change it in the script]
'''
