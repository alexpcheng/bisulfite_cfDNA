#!/bin/bash
#Script run to create multithreads.
#Assumes this script is run from project_name/Bin
SET=$1

#set amount of cores
NUM_PROC=51
mkdir -p ../V1/read_statistics
#build JOBLIST
cat ../lists/included_samples.txt | while read sample
do
	echo "bash read_mapping_statistics.sh $sample" >> prompt
done

#perl_fork_univ.pl prompt $NUM_PROC

cat prompt | while read line
do
	echo $line
done

echo -e "sample\ttotal_reads\thuman_aligned\tHQ_human_aligned\tmicrobe_mapped" > ../V1/read_statistics/BS_treated_samples.txt
cat ../V1/read_statistics/*.tmp >> ../V1/read_statistics/BS_treated_samples.txt
rm prompt
