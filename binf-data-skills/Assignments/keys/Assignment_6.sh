#!/bin/bash

#Assignment 6 key

#3.
wget https://gzahn.github.io/binf-data-skills/Data/Assignment_6_File_1.fasta.gz
gunzip Assignment_6_File_1.fasta.gz

#4. & 5.
paste <(zcat Assignment_6_File_1.fasta.gz | grep "^>" | sed 's/>//') <(zcat Assignment_6_File_1.fasta.gz | grep -v "^>") > Assignment_6_File_1.tsv

#6.
cat Assignment_6_File_1.tsv | cut -f 1 | sed 's/^OTU/>OTU/' > otulist
cat Assignment_6_File_1.tsv | cut -f 2 > seqlist
paste -d \\n otulist seqlist


#7.
wget https://gzahn.github.io/binf-data-skills/Data/Assignment_6_File_2.fasta.gz
gunzip Assignment_6_File_2.fasta.gz

#8.
cat Assignment_6_File_2.fasta | grep "Enterobacteriaceae" | sort -u | cut -d ";" -f 6

#9. & 10.
cat Assignment_6_File_2.fasta | grep "^>" | sort | grep "Enterobacteriaceae" > tmp
paste <(cat Assignment_6_File_2.fasta | grep "Enterobacteriaceae" | sort) <(while read line;do grep -c $line Assignment_6_File_2.fasta ;done <tmp) | sort -u > Enterobacteriaceae_Genus_Counts.tsv




