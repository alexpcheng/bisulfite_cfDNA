#!/bin/bash
# Title: get_HQ_lengths.sh
# Authors: Alexandre Pellan Cheng
# Brief description: Get the lengths of paired end sequencing data (mapQ >=40)

NUM_PROC=51
HQ_directory=../V1/aligned/HQ/
outdir=../V1/Lengths/
bam=_mapped_autosomal.hq.bam

cat ../lists/included_samples.txt | while read sample
do
	echo "samtools view $HQ_directory$sample$bam | awk '{if (\$9>0 && \$9<10000) print \$9}' > $outdir$sample.lengths" >> prompt
done

#perl_fork_univ.pl prompt $NUM_PROC #convenience function to multithread

cat prompt | while read line
do
	echo $line
done

rm prompt
