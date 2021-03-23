#!/bin/bash	


# download and move files to scratch storage space

wget -c -O /scratch/kingspeak/serial/u6033249/H7H77BGX2_1_fastq.gz "http://cbsuapps.tc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1353955503&refid=222144"

wget -q -c -O /scratch/kingspeak/serial/u6033249/HMHFMBGX2_1_fastq.gz "http://cbsuapps.biohpc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1651073382&refid=253718"

wget -q -c -O /scratch/kingspeak/serial/u6033249/HMJVMBGX2_1_fastq.gz "http://cbsuapps.biohpc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=227723316&refid=253717"

wget -q -c -O /scratch/kingspeak/serial/u6033249/HMG32BGX2_1_fastq.gz "http://cbsuapps.biohpc.cornell.edu/Sequencing/showseqfile.aspx?mode=http&cntrl=1334182982&refid=239744"


# make note of location
echo "/scratch/kingspeak/serial/u6033249/"
ls -ahl /scratch/kingspeak/serial/u6033249/


