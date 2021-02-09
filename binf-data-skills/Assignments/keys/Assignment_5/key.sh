#!/bin/bash
set -euo pipefail

# copy Case Study 1 summary file to this directory
cp ../../../Case_Studies/01_extract_ITS1_region/key/summary_output_per_file.tsv ./

# Proportion of Reads that passed ITSxpress
paste summary_output_per_file.tsv <(cat <(echo "PassingProportion") <(cat summary_output_per_file.tsv | tail -n +2 | awk '{print $3 / $2}')) > Passing_Proportions.tsv

# Find list of files that had fewer than 20% passing reads
cat Passing_Proportions.tsv | tail -n +2 | awk '$5 < 0.2 {print $1}' > Unsuccessful_ITSxpress_Files.txt

# Print total number of reads from those unsuccessful files (using Unsuccessful_Files...txt list)
while read line; do grep "^$line" summary_output_per_file.tsv; done < Unsuccessful_ITSxpress_Files.txt | awk 'BEGIN{s = 0}; {s += $2}; END{print s}'

# Using just awk on Passing_Porportions.tsv
cat Passing_Proportions.tsv | awk '$5 < 0.2' | awk 'BEGIN{s = 0}; {s += $2}; END{print s}'

# Find mean phred scores for each input file
paste <(ls -1 ../../../Case_Studies/01_extract_ITS1_region/key/DNA*fastq) <(for i in ../../../Case_Studies/01_extract_ITS1_region/key/DNA*fastq;do fastx_quality_stats -i $i | tail -n +2 | awk 'BEGIN{s = 0}; {s += $8}; END{print s/NR}';done) > Mean_phred_scores.tsv

# Add these results to full summary info table
paste Passing_Proportions.tsv <(cat <(echo "Mean_Quality") <(cat Mean_phred_scores.tsv | awk '{print $2}')) > Full_Summary_Output.tsv
