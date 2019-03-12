# bisulfite_cfDNA
Custom scripts used for data analysis of microbial and human derived cell-free DNA after bisulfite treatment published in [journal]

## Requirements
The following softwares need to be available, as they are called in certain parts of the pipeline
### Main pipeline
- snakemake
- Trimmomatic
- bwa-meth
- Samtools
- MethylDackel
- BEDTools
- BEDOPS
- methpipe

### Methylation references pipeline
- snakemake
- CrossMap
- Metilene
- BEDTools

### Genomic abundance of pathogens
- snakemake
- Samtools
- Trimmomatic
- FLASH
- SGREP
- bwa-meth
- bwa
- GRAMMy
## Folder structure
Most scripts are run relative to the current directory, so initialize your workspace similarly to the repo:
```
git clone https://github.com/alexpcheng/bisulfite_cfDNA/
```

## Main pipeline
```
# Index your reference genome using bwa-meth
bwameth.py index [REFERENCE]

# modify the config file with file names and locations of the previously mentioned software.

cd bisulfite_cfDNA
source activate [conda environment]
snakemake
```

## Methylation references pipeline
```
# modify the config file with file names and locations of the previously mentioned software.
cd bisulfite_cfDNA/Methylation_References_snaked
source activate [conda environment]
snakemake
```
## Genomic abundance of pathogens
```
# modify the config file with file names and locations of the previously mentioned software.
#make sure the BLASTDB variable is in your path and points to your NCBI blast database.
cd bisulfite_cfDNA/GRAMMy/
source activate [conda environment]
#prepare GRAMMy and BLAST databases
snakemake -s Snakefile.databases
# main pipeline
snakemake 
```
## Run scripts for specific tasks
```
cd bisulfite_cfDNA/Bin

# Donor fractions
bash donor_fractions_par.sh #make sure you have enough CPUs available, or change it in the script

# Get lengths profile of reads with high mapQ
bash get_HQ_lengths.sh

# Certain read mapping stats
bash read_mapping_statistics_par.sh
```

## Generate figures as found in [publication]
```
cd bisulfite_cfDNA/Bin
R..
```
