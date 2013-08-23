#!/usr/bin/perl -w
use strict;

# use this script below to classify them to normal pairs and single ends, 
# then you can shuffle pairs and merge the single ends for velvet input. 
# It is just available for Illumina reads. 


#http://biostar.stackexchange.com/questions/8747/illumina-pair-end-library-de-novo-assembly-broken-pairs

die "perl $0 <IN1:PE1.fastq> <IN2:PE2.fastq> <OUT1:out.PE1.fastq> <OUT2:out.PE2.fastq>\n" if (@ARGV!=4);

#@ILLUMINA-57021F:7:1:1020:11315#0/1
#@HWI-ST397_0000:6:1101:1455:2111#TTAGGC/1
my ($PE1,$PE2,$out1,$out2)=@ARGV;
my %hash;
open (IN,$PE1) || die $!;
while(<IN>)
{
    if (/^\@(\S+)\#\S+\/[12]/)
    {
        $hash{$1}=$_;
        $hash{$1}.=<IN>;
        $hash{$1}.=<IN>;
        $hash{$1}.=<IN>;
    }
}
close IN;

open (IN,$PE2) || die $!;
open (OUT1,">$out1") || die $!;
open (OUT2,">$out2") || die $!;
open (SE1,">$out1.single") || die $!;
open (SE2,">$out2.single") || die $!;
while(<IN>)
{
    if (/^\@(\S+)\#\S+\/[12]/)
    {
        if (exists $hash{$1})
        {
            print OUT1 $hash{$1};
            undef $hash{$1};
            delete $hash{$1};
            print OUT2 $_;
            $_=<IN>;
            print OUT2 $_;
            $_=<IN>;
            print OUT2 $_;
            $_=<IN>;
            print OUT2 $_;
        }
        else
        {
            print SE2 $_;
            $_=<IN>;
            print SE2 $_;
            $_=<IN>;
            print SE2 $_;
            $_=<IN>;
            print SE2 $_;
        }
    }
}
foreach my $key(keys %hash)
{
    if ((defined $hash{$key})&&($hash{$key} ne ""))
    {
        print SE1 $hash{$key};
    }
}

close IN;
close OUT1;
close OUT2;
close SE1;
close SE2;
exit;