#!/bin/sh -e

# Surya Saha
# Purpose:To write out paired fastq files with corresponding seqs in order
# Params: Fwd.fq Rev.fq Common_names

usage(){
	echo "Params: Fwd.fq Rev.fq"
	echo "Purpose:To write out paired fastq files with corresponding seqs in order"
}

if [ $# != 2 ]
then
	usage; exit 1;
fi

#get params
FWD=$1
REV=$2

#pre-processing
export LANG=C; export LC_ALL=C
fgrep --color=auto '@' $FWD  | sed 's,\/1,,'| sed 's,\@,,' | sort > 1.names.trimmed.sorted
fgrep --color=auto '@' $REV  | sed 's,\/2,,'| sed 's,\@,,' | sort > 2.names.trimmed.sorted
comm -12 1.names.trimmed.sorted 2.names.trimmed.sorted > commnames

CTR=0
while read line
do
	#write fwd read with 3 following lines, not found val = 1
	grep -A3 $line $FWD >> ordered.${FWD}
	if [ $? = "1" ] ; then
		printf "$line not found in $FWD\nTemp files not deleted.\n"; exit 1;
	fi
	
	#write rev read with 3 following lines
	grep -A3 $line $REV >> ordered.${REV}
	if [ $? = "1" ] ; then
		printf "$line not found in $REV\nTemp files not deleted.\n"; exit 1;
	fi
	
	CTR=$(($CTR+1));
done < commnames

#post-processing
rm -f 1.names.trimmed.sorted 2.names.trimmed.sorted commnames
echo "$CTR paired reads written.."

