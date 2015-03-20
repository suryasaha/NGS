#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Feb 8, 2011
# http://seqanswers.com/forums/showthread.php?t=10392
#!/usr/bin/env perl

# Program to compare two quality filtered fastq files from a PE run and output new files containing only shared sequences
# Input is quality filtered files for reads 1 and 2
# Output is two new fastq files, one for each read containing only shared reads
# Reads will be output in the same order in each file
# greigite 3.2011

if (@ARGV != 2) {    die "Usage: findUniqPairs.v2.pl <fastq for read 1> <fastq for read 2>\nNote: delete old shared files\n\n";}

print STDERR "Did you DELETE previously created .SHARED files??\nDoes NOT record UNPAIRED reads\n\n";

my (@readsone,@readstwo,@combined);
my $filenameone = $ARGV[0] . ".shared";
my $filenametwo = $ARGV[1] . ".shared";
#my $logfile = "unshared_reads.log";

open my $readoneout, ">>$filenameone";
open my $readtwoout, ">>$filenametwo";
#open my $logfh, ">>$logfile";

# build arrays of read names
my @readsonetmp = `grep "@" $ARGV[0]`;
my @readstwotmp = `grep "@" $ARGV[1]`;

@readsone = map{ chomp $_; $_ =~s/'//g; my ($trim) = $_ =~/(.*)\/[0-9]/; $trim;}@readsonetmp; 
@readstwo = map{ chomp $_; $_ =~s/'//g; my ($trim) = $_ =~/(.*)\/[0-9]/; $trim;}@readstwotmp; 

foreach my $readone(@readsone){
  my (@printonetrim,@printtwotrim);
  if (grep ($_ eq $readone, @readstwo)) {
     my $lineone = `grep -n $readone $ARGV[0]`;
     my ($linenum) = $lineone =~/([0-9]{1,}):.*/;
     $linenum = $linenum +3;
     my @printone = `head -$linenum $ARGV[0] | tail -4`;

     my $linetwo = `grep -n $readone $ARGV[1]`;
     ($linenum) = $linetwo =~/([0-9]{1,}):.*/;
     $linenum = $linenum +3;
     my @printtwo = `head -$linenum $ARGV[1] | tail -4`;

     print $readoneout @printone;
     print $readtwoout @printtwo;
  }
}

#=head1 NAME
#
#findUniqPairs.v2.pl - Remove duplicates reads from a file 
#
#=head1 SYNOPSIS
#
#  % findUniqPairs.v2.pl -fwd file -rev file -out file -combine file -informat Fasta/fastq -outformat Fasta/fastq
#  
#=head1 DESCRIPTION
#
#This script reads in paired read files and pulls out all unique pairs and writes out a file
#with combined reads. Orphans are written out to separate files.
#
#=head1 VERSION HISTORY
#
# Version 1: Standard
# Version 1.5: Can specify both input and output formats.
# Version 2: Completely new algo to reduce memory usage. 
#
#=head1 COMMAND-LINE OPTIONS
#
#Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
#are mandatory (see below).
#
#   --fwd       <>  File with fwd reads (required)
#   --rev       <>  File with rev reads (required)
#   --comb      <>  out file name (required)
#   --informat  <>  Fasta, genbank, EMBL, fastq (required) 
#   --outformat <>  Fasta, genbank, EMBL, fastq (required)
#      
#=head1 AUTHOR
#
#Surya Saha, ss2489@cornell.edu
#
#=cut
#
#my ($i,$j,$k,$fwd,$rev,$iformat,$oformat,$out,@temp,$ctr,%fwdctrs,%fwdseq,%fwdqual,
#%revseq,%revqual,%revctrs,);
#
#GetOptions (
#	'fwd=s' => \$fwd,
#	'rev=s' => \$rev,
#	'comb=s' => \$out,
#	'informat=s' => \$iformat,
#	'outformat=s' => \$oformat) or (system('pod2text',$0), exit 1);
#
## defaults and checks
#defined($fwd) or (system('pod2text',$0), exit 1);
#if (!(-e $fwd)){print STDERR "$fwd not found: $!\n"; exit 1;}
#defined($rev) or (system('pod2text',$0), exit 1);
#if (!(-e $rev)){print STDERR "$rev not found: $!\n"; exit 1;}
#if(($iformat ne 'Fasta') && ($iformat ne 'genbank') && ($iformat ne 'EMBL') && ($iformat ne 'fastq')){
#	system('pod2text',$0), exit 1;
#}
#if(($oformat ne 'Fasta') && ($oformat ne 'genbank') && ($oformat ne 'EMBL') && ($oformat ne 'fastq')){
#	system('pod2text',$0), exit 1;
#}
#
#my $infwd = Bio::SeqIO->new(-file=>$fwd, -format=>$iformat);
#my $inrev = Bio::SeqIO->new(-file=>$rev, -format=>$iformat);
#my $outfile = Bio::SeqIO->new(-file=>">$out", -format=>$oformat);
#my $fwdorph = Bio::SeqIO->new(-file=>">orphans.$fwd", -format=>$oformat);
#my $revorph = Bio::SeqIO->new(-file=>">orphans.$rev", -format=>$oformat);
#
##init ctrs
#$fwdctrs{'all'}=0; $fwdctrs{'uniq'}=0; $fwdctrs{'dup'}=0; $fwdctrs{'mate-p'}=0; $fwdctrs{'mate-ab'}=0;
#$revctrs{'all'}=0; $revctrs{'uniq'}=0; $revctrs{'dup'}=0; $revctrs{'mate-p'}=0; $revctrs{'mate-ab'}=0;
#
##reading in seqs
#while (my $obj = $infwd->next_seq()){
#	#$fwdseq{$obj->display_id()}=$obj;
#	$fwdseq{$obj->display_id()}='';
#	$fwdctrs{'all'}++;
#}
#$fwdctrs{'uniq'}=scalar (keys %fwdseq);
#$fwdctrs{'dup'}=$fwdctrs{'all'}-$fwdctrs{'uniq'};
#
#while (my $obj = $inrev->next_seq()){
##	$revseq{$obj->display_id()}=$obj;
#	$revseq{$obj->display_id()}='';
#	$revctrs{'all'}++;
#}
#$revctrs{'uniq'}=scalar (keys %revseq);
#$revctrs{'dup'}=$revctrs{'all'}-$revctrs{'uniq'};
#
##pairing fwds to revs
#while(($i,$j) = each %fwdseq){
#	$k=$i; $k=~ s/1$/2/;
#	if(exists $revseq{$k}){
#		if($oformat eq 'fastq'){
##			$fwdseq{$i}->quality_header=0;
#			$outfile->write_fastq($fwdseq{$i});
##			$outfile->write_qual($fwdseq{$i});
##			$revseq{$k}->quality_header(0);
#			$outfile->write_fastq($revseq{$k});
##			$outfile->write_qual($revseq{$k});
#		}
#		else{
#			$outfile->write_seq($fwdseq{$i});
#			$outfile->write_seq($revseq{$k});
#		}
#		$fwdctrs{'mate-p'}++;
#	}
#	else{
#		if($oformat eq 'fastq'){
##			$fwdseq{$i}->quality_header(0);
#			$fwdorph->write_fastq($fwdseq{$i});
##			$fwdorph->write_qual($fwdseq{$i});
#		}
#		else{
#			$fwdorph->write_seq($fwdseq{$i});
#		}
#		$fwdctrs{'mate-ab'}++;
#	}
#}
#
##counting for revs
#while(($i,$j) = each %revseq){
#	$k=$i; $k=~ s/2$/1/;
#	if(exists $fwdseq{$k}){
#		$revctrs{'mate-p'}++;
#	}
#	else{
#		if($oformat eq 'fastq'){
##			$revseq{$i}->quality_header(0);
#			$revorph->write_fastq($revseq{$i});
##			$revorph->write_qual($revseq{$i});
#		}
#		else{
#			$revorph->write_seq($revseq{$i});
#		}		
#		$revctrs{'mate-ab'}++;
#	}
#}
#
#if($fwdctrs{'mate-ab'}==0){ unlink "orphans.$fwd";}
#if($revctrs{'mate-ab'}==0){ unlink "orphans.$rev";}
#
##print counts
#print STDERR "\tSeqs\tUniq\tDup\tMate-present\tMate-absent\n";
#print STDERR "Fwd\t$fwdctrs{'all'}\t$fwdctrs{'uniq'}\t$fwdctrs{'dup'}\t$fwdctrs{'mate-p'}\t\t$fwdctrs{'mate-ab'}\n";
#print STDERR "Rev\t$revctrs{'all'}\t$revctrs{'uniq'}\t$revctrs{'dup'}\t$revctrs{'mate-p'}\t\t$revctrs{'mate-ab'}\n";
#my($user_t,$system_t,$cuser_t,$csystem_t); ($user_t,$system_t,$cuser_t,$csystem_t) = times;
#print STDERR "System time for process: ",sprintf("%.3f",$system_t/3600)," hrs\n"; print STDERR "User time for process: ",sprintf("%.3f",$user_t/3600)," hrs\n";
#
#exit;