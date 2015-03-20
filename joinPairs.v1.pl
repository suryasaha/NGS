#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha July 8, 2011

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Bio::SeqIO;
use lib "/home/surya/bin/modules";
use SS;

=head1 NAME

 joinPairs.v1.pl - Generate paired files from trimmed and fq-joined dataset 

=head1 SYNOPSIS

  % joinPairs.v3.pl -comm names --fwd file --rev file --joined file --ns 15 --informat Fasta/fastq
  
=head1 DESCRIPTION

 This script reads in paired read files that have been quality trimmed and overlapping reads joined so the order 
 has been disturbed. Produces a joined paired-read files and fwd/rev unpaired read files. Low memory requirement 
 but long runtime (4-6hrs for 30mill reads). 
 
 Accounting for out of order reads was unnecessary since fq-join files are in order. See bash commands
 
=head1 BASH COMMANDS
 export LANG=C; export LC_ALL=C
 fgrep --color=auto '@' s_6_1_sequence.txt.trim28.fastq  | sed 's,\/1,,'| sed 's,\@,,' | sort > 1.names.trimmed.sorted
 fgrep --color=auto '@' s_6_2_sequence.txt.trim28.fastq  | sed 's,\/2,,'| sed 's,\@,,' | sort > 2.names.trimmed.sorted
 comm -12 1.names.trimmed.sorted 2.names.trimmed.sorted > commnames

=head1 TODO
 Ver 2 Output in Fastq and Fasta using Bio::SeqIO
 
=head1 VERSION HISTORY

 Version 1: Shell commands + perl code. Order of paired reads is dependent on reads being in order in input files.
            Input: Common name file produced by comm, fwd/rev filtered paired end files
            Output: Joined ordered paired end files, fwd/rev unpaired files  

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --comm      <>  List of common read names (required)
   --fwd       <>  File with fwd reads (required)
   --rev       <>  File with rev reads (required)
   --joined    <>  Name for joined file (required)
   --ns        <>  Number of N's in join (default=15)
   --informat  <>  Fasta, genbank, EMBL, fastq (required) 
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$k,$fwd,$rev,$comm,$iformat,$n,$out,@temp,$rec,%cnames,@ctrs,%fData);

GetOptions (
	'comm=s'=> \$comm,
	'fwd=s' => \$fwd,
	'rev=s' => \$rev,
	'joined=s' => \$out,
	'ns:i' => \$n,
	'informat=s' => \$iformat) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($comm) or (system('pod2text',$0), exit 1);
if (!(-e $comm)){print STDERR "$comm not found: $!\n"; exit 1;}
defined($fwd) or (system('pod2text',$0), exit 1);
if (!(-e $fwd)){print STDERR "$fwd not found: $!\n"; exit 1;}
defined($rev) or (system('pod2text',$0), exit 1);
if (!(-e $rev)){print STDERR "$rev not found: $!\n"; exit 1;}
defined($out) or (system('pod2text',$0), exit 1);
$n ||= 15; 
if(($iformat ne 'Fasta') && ($iformat ne 'genbank') && ($iformat ne 'EMBL') && ($iformat ne 'fastq')){
	system('pod2text',$0), exit 1;
}

my $infwd = Bio::SeqIO->new(-file=>$fwd, -format=>$iformat);
my $inrev = Bio::SeqIO->new(-file=>$rev, -format=>$iformat);
unless(open(COMM,"<$comm")){print "not able to open $comm\n\n";exit 1;}
unless(open(OUTJ,">${out}")){print "not able to open joined.${out}\n\n";exit 1;}
unless(open(OUTFU,">unpaired.${fwd}.fas")){print "not able to open unpaired.${fwd}.fas\n\n";exit 1;}
unless(open(OUTRU,">unpaired.${rev}.fas")){print "not able to open unpaired.${rev}.fas\n\n";exit 1;}

#read in common names
while($rec=<COMM>){
	chomp $rec; $cnames{$rec}='';
}

#initing ctrs
for $i (0..3){$ctrs[$i]=0;} #fwd paired, fwd unpaired, rev paired, rev unpaired

print STDERR "NOTE: Output paired file order is dependent on order in input files.\n";
print STDERR "NOTE: Call fasta_formatter to format joined fasta file.\nStarting parse...\n";

#parsing fwd seqs
while (my $obj = $infwd->next_seq()){
	$i=$obj->display_id();
	$i=~ s/^\@//; $i=~ s/\/1$//;
	#print STDOUT "Read $i from fwd\n";#debug
	if(exists $cnames{$i}){
		$ctrs[0]++;
		$cnames{$i}=$obj->seq();#using cnames to remember fwd seq
		
#		if($oformat eq 'fastq'){
#			$outpfwd->write_fastq($obj);
#			#print STDOUT "Wrote $i to fwd paired\n";#debug
#		}
#		else{
#			$outpfwd->write_seq($obj);
#		}		
	}
	else{
		$ctrs[1]++;
		print OUTFU '>',$ctrs[1],"\n",$obj->seq(),"\n";
#		if($oformat eq 'fastq'){
#			$outufwd->write_fastq($obj);
#			#print STDOUT "\tWrote $i to fwd unpaired\n";#debug
#		}
#		else{
#			$outufwd->write_seq($obj);
#		}		
	}
}
print STDERR "Finished parsing Fwds..starting on Revs\n";
#parsing rev seqs
while (my $obj = $inrev->next_seq()){
	$i=$obj->display_id();
	$i=~ s/^\@//; $i=~ s/\/2$//;
	#print STDOUT "Read $i from rev\n";#debug
	if(exists $cnames{$i}){
		$ctrs[2]++;
		print OUTJ '>',$ctrs[2],"\n",$cnames{$i};
		for $j (1..$n){ print OUTJ 'N';}
		print OUTJ $obj->seq(),"\n";
		
#		if($oformat eq 'fastq'){
#			$outprev->write_fastq($obj);
#			#print STDOUT "Wrote $i to rev paired\n";#debug
#		}
#		else{
#			$outprev->write_seq($obj);
#		}		
	}
	else{
		$ctrs[3]++;
		print OUTRU '>',$ctrs[3],"\n",$obj->seq(),"\n";
#		if($oformat eq 'fastq'){
#			$outurev->write_fastq($obj);
#			#print STDOUT "\tWrote $i to rev unpaired\n";#debug
#		}
#		else{
#			$outurev->write_seq($obj);
#		}
	}
}
print STDERR "Finished parsing Revs\n\n";
close(OUTJ);
close(OUTFU);
close(OUTRU);

if($ctrs[0]!=$ctrs[2]){ print STDERR "Err in counts: fwd paired not equal rev paired\n\n"; exit 1;}

# print ctrs
print STDERR "Paired reads: $ctrs[0]\nUnpaired forwards: $ctrs[1]\nUnpaired reverses: $ctrs[3]\n\n";
# delete empty unpaired files
if($ctrs[1]==0){unlink "unpaired.${fwd}.fas";}
if($ctrs[3]==0){unlink "unpaired.${rev}.fas";}

&SS::mem_used();
&SS::runtime();
exit;