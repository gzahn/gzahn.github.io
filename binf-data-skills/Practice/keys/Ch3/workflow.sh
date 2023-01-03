#!/bin/bash

# get data
wget http://gzahn.github.io/binf-data-skills/Data/Chapter_3_Practice_Files.zip

# unzip
unzip Chapter_3_Practice_Files.zip

# 1.
paste Chapter_3_Practice_Files/A.txt Chapter_3_Practice_Files/B.txt > Chapter_3_Practice_Files/C.txt

# 2.
paste <(head -n 6 Chapter_3_Practice_Files/C.txt) <(grep -i "^E\|^B" Chapter_3_Practice_Files/words.txt)

# 3. Tricky!
cat Chapter_3_Practice_Files/story.txt | sed 's/\. /_/g' | cut -d "_" -f 1,2,3,4,5,6 | sed 's/_/\. /g'

# 4. Extra tricky!
paste <(head -6 Chapter_3_Practice_Files/A.txt) <(cat Chapter_3_Practice_Files/story.txt | sed 's/\. /_/g' | cut -d "_" -f 1,2,3,4,5,6 | sed 's/_/\n/g')

