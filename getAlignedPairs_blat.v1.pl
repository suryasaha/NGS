#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 8, 2011

use strict;
use warnings;
use Getopt::Long;
use POSIX;
eval {
	require Bio::SearchIO;
	require Bio::SeqIO;
};
use Bio::SearchIO;
use Bio::SeqIO;

=head1 NAME

getAlignedPairs_blat.v1.pl - Get the pairs which align to the reference genome with correct orientation
 

=head1 SYNOPSIS

  % getAlignedPairs.v1.pl -blatrep blat.out -pairedreads file.fasta -out file -outformat Fasta/fastq
  
=head1 DESCRIPTION

This script reads in blat output of paired read fasta file to reference genome and pulls out all 
pairs that align with correct orientation and spacing and writes out fwd and rev files with 
selected paired reads. Using hash of 2D arrays to store blat hits.

=head1 VERSION HISTORY
 Version 1: Standard

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --blatrep     <blat.out>    Blat report file (required)
   --matlen      <100>         Minimum match length(bp), if 100% match is reqd for 100 bp, then 100 (required)
   --pairedreads <file.fasta>  Input sequence file (required)
   --readlen     <100>         Length of input reads (required)
   --inslen      <int>         Insert length in bp (required)
   --stdev       <float>       Standard deviation of insert length in bp (required)
   --out         <file>        Output reads core file name (required)
   --format      <Fasta/fastq> Format of input/output files (required)

=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$k,$l,$m,$blatrep,$pairedreads,$rlen,$mlen,$inslen,$stdev,$format,$out,@temp,$flag,
@temp1,@temp2,%seqs,%hits,$fout,$rout,$minins,$maxins,$ctr);

GetOptions (
	'blatrep=s' => \$blatrep,
	'matlen=i' => \$mlen,
	'pairedreads=s' => \$pairedreads,
	'readlen=i' => \$rlen,
	'inslen=i' => \$inslen,
	'stdev=f' => \$stdev,
	'out=s' => \$out,
	'format=s' => \$format) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($blatrep) or (system('pod2text',$0), exit 1);
if (!(-e $blatrep)){print STDERR "$blatrep not found: $!\n"; exit 1;}
$out ||= "aligned.${pairedreads}";
if(($format ne 'Fasta') && ($format ne 'fastq')){
	system('pod2text',$0), exit 1;
}

$minins=$inslen-(2*$rlen)-$stdev;
$maxins=$inslen-(2*$rlen)+$stdev;

print STDERR "Using minimum match cutoff of $mlen bp ...\n";
print STDERR "Setting insert size range from $minins to $maxins bp \(negative value implies overlap of paired reads\)...\n";
if($format eq 'fastq'){print STDERR "Presuming $pairedreads is in $format format\n";}

#read in sequences
$i = Bio::SeqIO->new(-file=>"<$pairedreads", -format=>$format);
$k=0;
## POOR DESIGN SINCE MACHINE CHOKES ON RAM WHEN STORING AN OBJ FOR EACH SEQ ##
while ($j = $i->next_seq()){
	$seqs{$j->display_id()}=$j; $k++;
	if($k%100000 == 0){print STDERR "Read in $k sequences ...\n";}
}
print STDERR "\nRead in $k sequences ...\n";

#process report
#match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#100	0	0	0	0	0	0	0	+	HWI-EAS39X_10218_FC6210P7543951489#0/2	100	0	100	gi|190570478|ref|NC_010981.1|	1482455	1034057	1034157	1	100,	0,	1034057,
#100	0	0	0	0	0	0	0	+	HWI-EAS39X_10218_FC6210P7547041673#0/2	100	0	100	gi|190570478|ref|NC_010981.1|	1482455	896042	896142	1	100,	0,	896042,
#100	0	0	0	0	0	0	0	+	HWI-EAS39X_10218_FC6210P7518111916#0/2	100	0	100	gi|190570478|ref|NC_010981.1|	1482455	413640	413740	1	100,	0,	413640,

unless(open(IN,"<$blatrep")){print "not able to open $blatrep....\n\n";exit 1;}
for(1..5){$i=<IN>}
#using hash of arrays
while($i=<IN>){
	@temp=split("\t",$i);
	if($temp[0]>=$mlen){
		@temp1=();
		if(exists $hits{$temp[9]}){
			@temp1=@temp2=();
			$temp1[0]=$temp[8];#strand
			$temp1[1]=$temp[15];#start on ref genome
			$temp1[2]=$temp[16];#end on ref genome
			@temp2=$hits{$temp[9]};#get prev hit info and add to it
			push @temp2,[@temp1];	$hits{$temp[9]}=[@temp2];
		}
		else{
			@temp1=@temp2=();
			$temp1[0]=$temp[8];#strand
			$temp1[1]=$temp[15];#start on ref genome
			$temp1[2]=$temp[16];#end on ref genome
			push @temp2,[@temp1];	$hits{$temp[9]}=[@temp2];
		}
	}
}

$fout = Bio::SeqIO->new(-file=>">fwd.${out}", -format=>$format);
$rout = Bio::SeqIO->new(-file=>">rev.${out}", -format=>$format);

$ctr=0;
while(($i,$j) = each(%hits)){
	#get pair info, if eligible, then print reads
	#and delete pair in hash
	if($i=~ /\/1$/){
		$k=$i; $k=~ s/\/1$/\/2/;
		if(exists $hits{$k}){
			my(@fcoords,@rcoords);
			#jumping thru hoops to assign array refs
			@fcoords=@$j; $l=$hits{$k}; @rcoords=@$l;
			for $l (0..$#fcoords){
				for $m (0..$#rcoords){
					if($rcoords[$m][0] eq '+'){next;}#skip if mate hits + strand
					else{
						#if rev hit is within range
						if(($rcoords[$m][1]>=($fcoords[$l][2]+$minins)) && ($rcoords[$m][1]<=($fcoords[$l][2]+$maxins))){
							if($format eq 'fastq'){ $fout->write_fastq($seqs{$i}); $rout->write_fastq($seqs{$k});}
							else{$fout->write_seq($seqs{$i}); $rout->write_seq($seqs{$k});}
							#delete fwd and rev from hash
							delete $hits{$i}; delete $hits{$k};
							print STDERR '.'; last;
						}
					}
				}
			}
		}
	}
	elsif($i=~ /\/2$/){
		$k=$i; $k=~ s/\/2$/\/1/;
		if(exists $hits{$k}){
			my(@fcoords,@rcoords);
			#jumping thru hoops to assign array refs
			@rcoords=@$j; $l=$hits{$k}; @fcoords=@$l;
			for $l (0..$#fcoords){
				for $m (0..$#rcoords){
					if($rcoords[$m][0] eq '+'){next;}#skip if mate hits + strand
					else{
						#if rev hit is within range
						if(($rcoords[$m][1]>=($fcoords[$l][2]+$minins)) && ($rcoords[$m][1]<=($fcoords[$l][2]+$maxins))){
							if($format eq 'fastq'){ $fout->write_fastq($seqs{$i}); $rout->write_fastq($seqs{$k});}
							else{$fout->write_seq($seqs{$i}); $rout->write_seq($seqs{$k});}
							#delete fwd and rev from hash
							delete $hits{$i}; delete $hits{$k};
							print STDERR '.'; last;
						}
					}
				}
			}
		}
	}
	$ctr++;
}

print STDERR "\nWrote out in $ctr paired read(s)\n";
exit;