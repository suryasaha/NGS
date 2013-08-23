#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 5, 2011

use strict;
use warnings;
use Getopt::Long;
use POSIX;
eval {
	require Bio::SeqIO;
};
use Bio::SeqIO;

=head1 NAME

removeDups.v1.pl - Remove duplicates reads from a file 

=head1 SYNOPSIS

  % removeDups.v1.pl 
  
=head1 DESCRIPTION

This script reads in a sequence file and pulls out the duplicates. Duplicates are 
determined from name or sequence similarity. Originally designed for Justins reads

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --reads     <>  File with reads (required)
   --informat  <>  Fasta, genbank, EMBL, fastq (required)
   --compare   <>  seq, name (required)
   --outformat <>  Fasta, genbank, EMBL, fastq 
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$reads,$informat,$compare,$outformat,@temp,$ctr);

GetOptions (
	'reads=s' => \$reads,
	'informat=s' => \$informat,
	'compare=s' => \$compare,
	'outformat:s' => \$outformat) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($reads) or (system('pod2text',$0), exit 1);
if (!(-e $reads)){print STDERR "$reads not found: $!\n"; exit 1;}
if(($informat ne 'Fasta') && ($informat ne 'genbank') && ($informat ne 'EMBL') && ($informat ne 'fastq')){
	system('pod2text',$0), exit 1;
}
if(($compare ne 'seq') && ($compare ne 'name')){ system('pod2text',$0), exit 1;}
$outformat ||=$informat;
if(($outformat ne 'Fasta') && ($outformat ne 'genbank') && ($outformat ne 'EMBL')&& ($informat ne 'fastq')){
	system('pod2text',$0), exit 1;
}

#my $in = Bio::SeqIO->new(-format => 'fastq-illumina', -file => $fwd);
my $in = Bio::SeqIO->new(-file=>$reads, -format=>$informat);
my $out = Bio::SeqIO->new(-file=>">uniq.$reads", -format=>$outformat);
my $dup = Bio::SeqIO->new(-file=>">dup.$reads", -format=>$outformat);

my %matching_hash = ();
$i=0;
while (my $obj = $in->next_seq()){
	if ($compare eq 'seq'){
		if(!$matching_hash{$obj->seq()}){
			if ($outformat eq 'fastq'){$out->write_fastq($obj);}
			else{$out->write_seq($obj);}
			#$final_hash{$obj->display_id} = $obj->seq;
			$matching_hash{$obj->seq()} = 1;
		}
		else{
			if ($outformat eq 'fastq'){$dup->write_fastq($obj);}
			else{$dup->write_seq($obj);}
		}
	}
	elsif($compare eq 'name'){
		if(!$matching_hash{$obj->display_id()}){
			if ($outformat eq 'fastq'){$out->write_fastq($obj);}
			else{$out->write_seq($obj);}
			#$final_hash{$obj->display_id} = $obj->seq;
			$matching_hash{$obj->display_id()} = 1;
		}
		else{
			if ($outformat eq 'fastq'){$dup->write_fastq($obj);}
			else{$dup->write_seq($obj);}
		}
	}
	$i++;
}

$j=scalar(values(%matching_hash));
print STDERR "Original nof seq: ",$i," \($reads\)\nUnique seqs: ",$j," \(uniq\.$reads\)\nDuplicate seqs: ",
	$i-$j," \(dup\.$reads\)\n";

exit;