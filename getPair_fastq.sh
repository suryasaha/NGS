#!/bin/bash
# Surya Saha
# SOL Genomics @BTI / Plant Path @Cornell
# Purpose: Fish out pair given fwd/rev of a paired end sequence

set -o nounset
set -o errexit

usage(){
	echo "Usage: $0 bait.fastq data.fastq"
	exit 1
}

log() { # classic logger
  local prefix="[$(date +%Y/%m/%d\ %H:%M:%S)]: "
  echo "${prefix} $@" >&2
}

WDIR=`pwd`
BAIT=$1 #fastq, bait pair
DATA=$2 #fastq, full set

grep '^@' "$BAIT" | sed 's,\/[12]$,,' > ${BAIT}.names
grep -A3 -f ${BAIT}.names "$DATA" > pairs.${BAIT}
rm -f ${BAIT}.names

log

 

