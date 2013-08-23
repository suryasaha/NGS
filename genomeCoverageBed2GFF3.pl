#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Jan 26, 2012

use strict;
use warnings;
use Getopt::Long;

=head1 NAME

 genomeCoverageBed2GFF3.pl - Convert genomeCoverageBed2GFF3 output to GFF 

=head1 SYNOPSIS

  % genomeCoverageBed2GFF3.pl --infile in.dz.out
  
=head1 DESCRIPTION

 Writes out covered regions with color code for coverage
 "genomeCoverageBed -ibam 10read.e2e.sorted.bam -g CP000857.1.genome -dz > 10.dz.out"

=head1 VERSION HISTORY

=head1 COMMAND-LINE OPTIONS

 Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --infile <.dz.out>    genomeCoverageBed output (required)
      
=head1 AUTHOR

 Surya Saha, ss2489@cornell.edu

=cut

my ($i,$rec,$in,$j,@temp,$name,$cov,$st,$end,$flag);

GetOptions (
	'infile=s' => \$in) or (system('pod2text',$0), exit 1);
defined($in) or (system('pod2text',$0), exit 1);

unless(open(IN,"<$in")){print "not able to open $in\n\n";exit 1;}
unless(open(OUT,">$in\.gff")){print "not able to open $in\.gff\n\n";exit 1;}
print OUT "\#\#gff-version 3\n";
#gi|224466365|gb|CP000857.1|	2796239	9
#gi|224466365|gb|CP000857.1|	2796240	9
#gi|224466365|gb|CP000857.1|	2796241	10
#gi|224466365|gb|CP000857.1|	2796242	10
#gi|224466365|gb|CP000857.1|	2796243	10
#gi|224466365|gb|CP000857.1|	2796244	10
#gi|224466365|gb|CP000857.1|	2796245	10
#gi|224466365|gb|CP000857.1|	2796246	11
#gi|224466365|gb|CP000857.1|	2796247	11
$flag=0; 
while($rec=<IN>){
	chomp $rec; @temp=split("\t",$rec);
	if($flag==0){
		$st=$end=$temp[1]+1;#since coords are 0 based
		$cov=$temp[2];
		$flag=1;
		$name=$temp[0];
	}
	elsif(($cov==$temp[2]) && (($end+1)==($temp[1]+1)))#same cov and next base
	{
		$end++; next;
	}
	elsif(($cov!=$temp[2]) || (($end+1)!=($temp[1]+1)))#diff cov or non-contigous base
	{
		print OUT "$name\tBEDTools\tmisc_feature\t$st\t$end\t$cov\t\.\t\.\tcolour\=";
		#20X red, 10X light red, 5X blue, 3X light blue, 1X light grey
		if($cov >= 20 ){print OUT "2\;\n";}
		elsif($cov >= 10 ){print OUT "16\;\n";}
		elsif($cov >= 5 ){print OUT "4\;\n";}
		elsif($cov >= 3 ){print OUT "9\;\n";}
		elsif($cov >= 1 ){print OUT "13\;\n";}
		#for next interval
		$st=$end=$temp[1]+1;#since coords are 0 based
		$cov=$temp[2];
	}
}
#for last record
print OUT "$name\tBEDTools\tmisc_feature\t$st\t$end\t$cov\t\.\t\.\tcolour\=";
		#20X red, 10X light red, 5X blue, 3X light blue, 1X light grey
		if($cov >= 20 ){print OUT "2\;\n";}
		elsif($cov >= 10 ){print OUT "16\;\n";}
		elsif($cov >= 5 ){print OUT "4\;\n";}
		elsif($cov >= 3 ){print OUT "9\;\n";}
		elsif($cov >= 1 ){print OUT "13\;\n";}
close(IN); close(OUT);
exit;
