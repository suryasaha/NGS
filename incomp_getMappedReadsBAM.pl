#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 6, 2012

use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Sam;
 
=head1 NAME

 getMappedReadsBAM.pl - Get uniquely and commonly mapped reads for a pair of BAM files 

=head1 SYNOPSIS

  % getMappedReadsBAM.pl --1inbam 1in.bam --2inbam 2in.bam
  
=head1 DESCRIPTION

 Writes out list of read names in 3 files

=head1 VERSION HISTORY

=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --1inbam <.bam>    1st BAM file (required)
   --2inbam <.bam>    2nd BAM file (required)
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

my ($i,$rec,$in1,$in2,$j,@temp,$flag,$cov,$st,$end,$len);

GetOptions (
	'1inbam=s' => \$in1,
	'2inbam=s' => \$in2) or (system('pod2text',$0), exit 1);
defined($in1) or (system('pod2text',$0), exit 1);
defined($in2) or (system('pod2text',$0), exit 1);

#unless(open(IN,"<$in")){print "not able to open $in\n\n";exit 1;}
unless(open(O1,">$in1\.uniq\.names")){print "not able to open $in1\.uniq\.names\n\n";exit 1;}
unless(open(O2,">$in2\.uniq\.names")){print "not able to open $in2\.uniq\.names\n\n";exit 1;}
unless(open(OC,">$in1\.$in2\.common\.names")){print "not able to open $in1\.$in2\.common\.names\n\n";exit 1;}


# high level API
#my $bam1 = Bio::DB::Sam->new(-bam  =>"data/ex1.bam", -fasta=>"data/ex1.fa",);
my $bam1 = Bio::DB::Sam->new(-bam  =>"data/ex1.bam",);



close(O1); close(O2); close(OC);


exit;