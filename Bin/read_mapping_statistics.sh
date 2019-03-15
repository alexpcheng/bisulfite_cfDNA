#!/bin/bash
# Title: read_mapping_statistics.sh
# Authors: Alexandre Pellan Cheng
# Brief description: Collect a few read statistics

sample=$1

#total number of reads
all_reads=$(zcat ../Data/$sample"_R1.fastq.gz" | wc -l)
number_of_reads=$(( all_reads / 4 ))

#number of reads that mapped back to hg19
bamfile=../V1/aligned/$sample.bam
raw_amount=$(samtools view -f 0x2 -c $bamfile)
r_a=$(( raw_amount / 2 ))
HQ_amount=$(samtools view -f 0x2 -q 40 -c $bamfile)
h_a=$(( HQ_amount / 2 ))

#number of reads that mapped back to a bacterial genome
grammy=../V1/infections/$sample.grammy.tab
number_of_microbial_reads=$(Rscript total_blast_hits.R $grammy)

echo -e "$sample\t$number_of_reads\t$r_a\t$h_a\t$number_of_microbial_reads" > ../V1/read_statistics/$sample.tmp
