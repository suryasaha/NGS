#!/bin/sh

# Surya Saha
# SRC http://seqanswers.com/forums/showthread.php?t=10392&goto=newpost
# Purpose: Merge fwd and revs into a single file.


#!/bin/bash
echo "Step #1 running...Estracting read names"
echo "Estracting forward read names"
grep "@HWI" basename.f.fastq > basename.list
echo "Estracting reverse read names"
grep "@HWI" basename.r.fastq >> basename.list

echo "Step #3 running...Estracting paired reads mames"
sed s/"\/[0-9]"//g basename.list | sort | uniq -d > basename.list.duplicated

echo "Step #3 running...Merging forward and reverse paired reads"
array=(`cat basename.list.duplicated`)
len=${#array[*]}

i=0
while [ $i -lt $len ]; do
       grep -w -A 3 "${array[$i]}" basename.f.fastq >> basename.merged.fastq
       grep -w -A 3 "${array[$i]}" basename.r.fastq >> basename.merged.fastq
       let i++
done
