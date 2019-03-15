#!/bin/bash
#Script run to create multithreads.
#Assumes this script is run from project_name/Bin
SET=$1

#set amount of cores
NUM_PROC=20

#Create list of patients who are sex mismatched
cat "../lists/recipient_donor_sexes.txt" | tail -n +2 | while read line
do
	echo "$line" | awk -F'\t' '{if ($3 != $2) print $0'} >> "../lists/only_sexmm.list"
done

#build JOBLIST
cat ../lists/only_sexmm.list | while read sample
do
	echo "bash donor_fractions.sh $sample" >> prompt 
done

#perl_fork_univ.pl prompt $NUM_PROC convenience function to multithread

cat prompt | while read line
do
	echo $line
done

Rscript donor_fractions_make_list.R
rm prompt
rm ../lists/only_sexmm.list
