#!/bin/sh

# Surya Saha
# BTI/PPath@Cornell
# Purpose: Write out a GFF3 file for a BED file supplied as parameter

#SL2.50ch03	18370591	18592711	Contig89
#SL2.50ch03	SL2.50_assembly	region	6521993	6866487	.	+	.	ID=1

printf "##gff-version 3\n"
awk -F"\t" 'BEGIN{OFS="\t"} {print $1,"source","region",$2+1,$3,".","+",".",$4}' "$1"



