#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Mar 15, 2011

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Bio::SeqIO;

=head1 NAME

 combinePairs.v1.pl - Combine 2 paired read files that are in same order 

=head1 SYNOPSIS

  % combinePairs.v1.pl -fwd file -rev file -out file -combine file -informat Fasta/fastq -outformat Fasta/fastq
  
=head1 DESCRIPTION

 This script reads in paired read files that are in the same order and writes out a file
 with combined reads. Stops if reads found out of order or in number of reads in both files
 is not equal. Designed for Illumina reads. Same functionality in ShuffleSequences.pl that
 comes with Velvet :-) but thats only for Fasta files

=head1 VERSION HISTORY

 Version 1: Standard

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --fwd       <>  File with fwd reads (required)
   --rev       <>  File with rev reads (required)
   --comb      <>  out file name (required)
   --informat  <>  Fasta, genbank, EMBL, fastq (required) 
   --outformat <>  Fasta, genbank, EMBL, fastq (required)
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$k,$fwd,$rev,$iformat,$oformat,$out,@temp,$ctr,$flag,$fobj,$robj);

GetOptions (
	'fwd=s' => \$fwd,
	'rev=s' => \$rev,
	'comb=s' => \$out,
	'informat=s' => \$iformat,
	'outformat=s' => \$oformat) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($fwd) or (system('pod2text',$0), exit 1);
if (!(-e $fwd)){print STDERR "$fwd not found: $!\n"; exit 1;}
defined($rev) or (system('pod2text',$0), exit 1);
if (!(-e $rev)){print STDERR "$rev not found: $!\n"; exit 1;}
if(($iformat ne 'Fasta') && ($iformat ne 'genbank') && ($iformat ne 'EMBL') && ($iformat ne 'fastq')){
	system('pod2text',$0), exit 1;
}
if(($oformat ne 'Fasta') && ($oformat ne 'genbank') && ($oformat ne 'EMBL') && ($oformat ne 'fastq')){
	system('pod2text',$0), exit 1;
}
if(($oformat eq 'fastq') && ($iformat ne 'fastq')){ print STDERR "Input format needs to be fastq if output format is fastq\n"; exit 1;}

my $infwd = Bio::SeqIO->new(-file=>$fwd, -format=>$iformat);
my $inrev = Bio::SeqIO->new(-file=>$rev, -format=>$iformat);
my $outfile = Bio::SeqIO->new(-file=>">$out", -format=>$oformat);


$ctr=0; $flag=1;
while($flag){
	$fobj = $infwd->next_seq();
	$robj = $inrev->next_seq();
	if((!$fobj) && (!$robj)){ print STDERR "Completed!\n"; last;}
	elsif(($fobj) && (!$robj)){print STDERR "ERR: $rev has less reads than $fwd..exciting\n"; last;}
	elsif((!$fobj) && ($robj)){print STDERR "ERR: $fwd has less reads than $rev..exciting\n"; last;}
	elsif(($fobj) && ($robj)){
		$i=$fobj->display_id(); $j=$robj->display_id();
		$i =~ s/1$//; $j =~ s/2$//; #remove trailing identifiers
		if($i ne $j){print STDERR $fobj->display_id(),' and ',$robj->display_id(),' are not paired..exciting',"\n"; last;}
		else{
			if($oformat eq 'fastq'){
				$outfile->write_fastq($fobj);
				$outfile->write_fastq($robj);
			}
			else{
				$outfile->write_seq($fobj);
				$outfile->write_seq($robj);
			}
		}
		$ctr++;
	}
}

print STDERR "Combined $ctr read pairs\n";
my($user_t,$system_t,$cuser_t,$csystem_t); ($user_t,$system_t,$cuser_t,$csystem_t) = times;
print STDERR "System time for process: ",sprintf("%.3f",$system_t/3600)," hrs\n"; print STDERR "User time for process: ",sprintf("%.3f",$user_t/3600)," hrs\n";

exit;