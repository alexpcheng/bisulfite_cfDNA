# bisulfite_cfDNA
Custom scripts used for data analysis of microbial and human derived cell-free DNA after bisulfite treatment published in [...]
Major data analysis pipeline is written using the snakemake workflow (https://snakemake.readthedocs.io/en/stable/).

## Required pipelines
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
- CrossMap
- Metilene
- BEDTools

### Genomic abundance of pathogens
- Samtools
- Trimmomatic
- FLASH
- SGREP
- bwa-meth
- bwa

## Creating databases, indexing references...

## Folder structure
Most scripts are run relative to the current directory, so initialize your workspace similarly to the repo:
```
git clone https://github.com/alexpcheng/bisulfite_cfDNA/
```
## Main pipeline
```
cd bisulfite_cfDNA
source activate [conda environment]
snakemake
```

## Methylation references pipeline
```
cd bisulfite_cfDNA/Methylation_References_snaked
source activate [conda environment]
snakemake
```
## Genomic abundance of pathogens
```
#make sure the BLASTDB variable is in your path and contains the following folder: 
# ...
cd bisulfite_cfDNA/GRAMMy/
snakemake 
```
## Run scripts for specific tasks
```
cd bisulfite_cfDNA/Bin

# Donor fractions
bash donor_fractions_par.sh #make sure you have enough CPUs available, or change it in the script

# Reads per sample
bash fastq_reads_par.sh

# Get lengths profile of reads with high mapQ
bash get_HQ_lengths.sh

# Certain read mapping stats
bash read_mapping_statistics_par.sh
```

## Generate figures as found in [publication]
```
```
