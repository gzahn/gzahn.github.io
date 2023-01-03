#!/bin/bash
set -e
set -u
set -o pipefail


find ../.. -name "Assignment_3_File_*fasta" > tmp


while read line;
do echo $line;
sleep 5;
echo "^^^That's the filename I'm looking at right now...";
sleep 5;
echo "Here are the first 2 lines of that file...";
sleep 5;
head -2 $line;
sleep 5;
echo "There are $(wc -l $line | cut -d " " -f 1) total lines in that file";
sleep 5;
echo "I'm about to show you each line in the file that contains GAACC..."
sleep 5;
grep "GAACC" $line;
sleep 5;
echo "Pretty neat, right!?"
sleep 2;
echo "Okay, now I'm going to move on to the next file..."
sleep 5;

done < tmp
rm -f tmp
echo "Oops, it looks like I'm out of files. That's the whole program." 
echo "I hope you enjoyed it. Please run it again any time you like!"
