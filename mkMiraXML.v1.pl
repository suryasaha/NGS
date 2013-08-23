#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Apr 29, 2011

use strict;
use warnings;
use Getopt::Long;
use POSIX;
use Bio::SeqIO;

=head1 NAME

 mkMiraXML.v1.pl - Create a NCBI trace XML file for a PE or MP solexa library for Mira 

=head1 SYNOPSIS

  % mkMiraXML.v1.pl -reads file -insert <int> -stdev <float> -xml out.xml -format Fasta/fastq -header_footer 0/1
  
=head1 DESCRIPTION

 This script reads in combined paired read files and writes out a NCBI TRACE-INFO format. 
 Designed for Illumina reads. 

=head1 VERSION HISTORY

 Version 1: Standard

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --reads         <>      File with combined reads (required)
   --insert        <>      Insert length (int,required)
   --stdev         <>      Standard dev  (int,required)
   --xml           <>      NCBI TRACE-INFO XML file
   --format        <>      Fasta, genbank, EMBL, fastq (required) 
   --header_footer <0/1>   Add header/footer info? (0 or 1, required)
      
=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$k,$reads,$xml,$format,$ins,$stdev,$hf,@temp,$ctr,$flag,$fobj,$robj);

GetOptions (
	'reads=s' => \$reads,
	'insert=i' => \$ins,
	'stdev=i' => \$stdev,
	'xml=s' => \$xml,
	'format=s' => \$format,
	'header_footer=i' => \$hf) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($reads) or (system('pod2text',$0), exit 1);
if (!(-e $reads)){print STDERR "$reads not found: $!\n"; exit 1;}
$xml ||= "$reads\.xml";
if(($format ne 'Fasta') && ($format ne 'genbank') && ($format ne 'EMBL') && ($format ne 'fastq')){
	system('pod2text',$0), exit 1;
}

unless(open(XML,">$xml")){print "not able to open $xml\n\n";exit 1;}

if($hf){print XML '<?xml version="1.0"?>',"\n",'<trace_volume>',"\n";}

my $in = Bio::SeqIO->new(-file=>$reads, -format=>$format);

while($i=$in->next_seq()){
	print XML '<trace>',"\n",'<trace_name>';
	print XML $i->display_id();
	print XML '</trace_name>',"\n",'<insert_size>',$ins,'</insert_size>',"\n";
	print XML '<insert_stdev>',$stdev,'</insert_stdev>',"\n",'</trace>',"\n";
}

if($hf){print XML '</trace_volume>',"\n";}
close(XML);
exit;