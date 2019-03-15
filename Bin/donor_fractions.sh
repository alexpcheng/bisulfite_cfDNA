#!/bin/bash
# Title: Donor fraction measurements in sex mismatched patients
# Authors: Philip Burnham, Alexandre Pellan Cheng
# Brief description: Called by donor_fractions_par.sh, makes donor fraction measurements.

# Assumes that the file is run from project_name/Bin

sample=$1
path="../V1/donor_fractions"
mkdir -p $path

if [ ! -f ../V2/donor_fractions/$sample.M.autoY.sexmm ]
then
  samtools view -b -q 40 ../V1/aligned/$sample"_mapped_all_chr.bam" > $path/$sample"_all_chr_hq.sorted.bam"
  samtools index $path/$sample"_all_chr_hq.sorted.bam"
  ../software/HMMcopy/bin/readCounter -w 500 \
	  		-c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
			$path/$sample"_all_chr_hq.sorted.bam" > \
  			$path/$sample".auto.readcounts.wig"
  Rscript donor_fractions_calculation.R $sample
fi
