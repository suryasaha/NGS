#!/usr/bin/perl -w
# PPath
# Surya Saha DATE??
# Purpose: 

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Bio::SeqIO::fastq;

=head1 NAME

removeDups.v1.pl - Remove duplicates reads from a file 

=head1 SYNOPSIS

  % removeDups.v1.pl 
  
=head1 DESCRIPTION

This script reads in 2 read files and pulls out the duplicates. Originally designed for Justins reads

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --reads     <>  File with reads (required)
   --informat  <>  Fasta, genbank, EMBL (required)
   --outformat <>  Fasta, genbank, EMBL (required)
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$reads,$informat,$outformat,@temp,$ctr);

GetOptions (
	'reads=s' => \$reads,
	'informat=s' => \$informat,
	'outformat:s' => \$outformat) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($reads) or (system('pod2text',$0), exit 1);
if (!(-e $reads)){print STDERR "$reads not found: $!\n"; exit 1;}
$outformat ||=$informat;

#my $in = Bio::SeqIO->new(-format => 'fastq-illumina', -file => $fwd);
my $in = Bio::SeqIO->new(-file=>$reads, -format=>$informat);
my $out = Bio::SeqIO->new(-file=>">uniq.$reads", -format=>$outformat);
my $dup = Bio::SeqIO->new(-file=>">dup.$reads", -format=>$outformat);

my %matching_hash = ();
$i=0;
while (my $obj = $in->next_seq()){
	if(!$matching_hash{$obj->seq()}){
		#$out->write_seq($obj);
		$out->write_fastq($obj);
		#$final_hash{$obj->display_id} = $obj->seq;
		$matching_hash{$obj->seq()} = 1;
	}
	else{
		#$dup->write_seq($obj);
		$dup->write_fastq($obj);
	}
	$i++;
}

$j=scalar(values(%matching_hash));
print $i-$j," duplicate(s) found\n";

exit;
