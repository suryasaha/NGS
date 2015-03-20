#!/usr/bin/perl -w
# PPath@Cornell
# Surya Saha Apr 11, 2011

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

getAlignedPairs_blastn.v1.5.pl - Get the pairs which align to the reference genome with correct orientation
 
=head1 SYNOPSIS

getAlignedPairs.v1.5.pl -blastnrep blastn.out -pairedreads file -matepair 0/1 -out file -outformat Fasta/fastq
  
=head1 DESCRIPTION

This script reads in blast output of paired read fasta file to reference genomes and pulls out all 
pairs that align with correct orientation and spacing. It then goes through the paired reads file 
and pulls out all qualifying read pairs and writes out selected paired reads. Also writes out a XLS 
file of mapping info. Using hash of 2D arrays to store blast hits. 

If matepair is true, use REVCOMPED files created using fastx_reverse_complement. It then also looks for 
pairs that map to ref genome as outies. These were the PE contaminants in the MP lib. These are written 
out with the PEcontam prefix.

=head1 VERSION HISTORY
 
Version 1: Only pulls out reads that match to wolbachia genomes. Handles only innies, i.e. short insert PE 
libs or revcomped MP libs without any PE contamination.
Version 1.5: Can deal with revcomped MP libs with PE contamination. Fixed a bug since e>s for all hits
irrespective of strand but I had presumed s>e for -1 strand.

=head1 TODO

1. Figure out how to use outlier and other endosymbiont results
2. How to screen against host genome?

=head1 COMMAND-LINE OPTIONS

Command-line options can be abbreviated to single-letter options, e.g. -f instead of --file. Some options
are mandatory (see below).

   --blastnrep   <blastn.out>  Blastn report file (required)
   --matlen      <100>         Minimum match length(bp), if identical match is reqd for 100 bp, then 100 (required)
   --pairedreads <file>        Input sequence file (required)
   --matepair    <0/1>         Input generated from Mate Pair lib? Must be REVCOMPED (required)
   --readlen     <100>         Length of input reads (required)
   --inslen      <int>         Insert length in bp (includes read length, required)
   --stdev       <int>         Standard deviation of insert length in bp (required)
   --out         <file>        Output reads core file name (required)
   --format      <Fasta/fastq> Format of input/output files (required)

=head1 AUTHOR

Surya Saha, ss2489@cornell.edu

=cut

my ($i,$j,$k,$l,$m,$blastnrep,$pairedreads,$rlen,$mlen,$inslen,$stdev,$format,$out,@temp,$flag,
$mp,@temp1,$minins,$maxins,$ctr,%genomes,$debug);

#remember cmd line..KLUDGE!!
$j=''; foreach $i (@ARGV){ $j.="$i ";}

GetOptions (
	'blastnrep=s' => \$blastnrep,
	'matlen=i' => \$mlen,
	'pairedreads=s' => \$pairedreads,
	'matepair=i' => \$mp,
	'readlen=i' => \$rlen,
	'inslen=i' => \$inslen,
	'stdev=f' => \$stdev,
	'out=s' => \$out,
	'format=s' => \$format,
	'debug:i' => \$debug,) or (system('pod2text',$0), exit 1);

# defaults and checks
defined($blastnrep) or (system('pod2text',$0), exit 1);
if (!(-e $blastnrep)){print STDERR "$blastnrep not found: $!\n"; exit 1;}
$out ||= "aligned.${pairedreads}";
if(($format ne 'Fasta') && ($format ne 'fastq')){
	system('pod2text',$0), exit 1;
}
unless(open(XLS,">mappingInfo.${out}.xls")){print "not able to open mappingInfo.${out}.xls\n\n";exit 1;}
print XLS "\t\tMAPPING REPORT\n\n\nParams: $j\n"; 
print XLS "\nCombined reads:\t$out\nBlastn:\t$blastnrep\nFwd reads:\tfwd.wol.${out}\nRev reads:\trev.wol.${out}\n";
if($mp){print XLS "Fwd PEcontaminant reads:\tfwd.PEcontam.wol.${out}\nRev PEcontaminant reads:\trev.PEcontam.wol.${out}\n";}

$minins=$inslen-(2*$rlen)-$stdev;#now insert size does not include read lengths
$maxins=$inslen-(2*$rlen)+$stdev;

print STDERR "Minimum match cutoff: $mlen bp ...\n";
print STDERR "Insert size range: $minins to $maxins bp \(exluding read lengths, negative value implies overlap of paired reads\)...\n";
if($format eq 'fastq'){print STDERR "Presuming $pairedreads is in $format format\n";}
if($mp){ print STDERR "Mate pair lib must have been REVCOMPED before blastn!!\n\n";}

#reference and outlier genome list
%genomes = (
'gi|50083297|ref|NC_005966.1| Acinetobacter sp. ADP1 chromosome, complete genome' => 'dce',
'gi|213155358|ref|NC_011585.1| Acinetobacter baumannii AB0057 plasmid pAB0057, complete sequence' => 'dce',
'gi|213155370|ref|NC_011586.1| Acinetobacter baumannii AB0057, complete genome' => 'dce',
'gi|215481761|ref|NC_011595.1| Acinetobacter baumannii AB307-0294, complete genome' => 'dce',
'gi|184159988|ref|NC_010605.1| Acinetobacter baumannii ACICU plasmid pACICU1, complete sequence' => 'dce',
'gi|184160017|ref|NC_010606.1| Acinetobacter baumannii ACICU plasmid pACICU2, complete sequence' => 'dce',
'gi|184156320|ref|NC_010611.1| Acinetobacter baumannii ACICU, complete genome' => 'dce',
'gi|126640097|ref|NC_009083.1| Acinetobacter baumannii ATCC 17978 plasmid pAB1, complete sequence' => 'dce',
'gi|126640109|ref|NC_009084.1| Acinetobacter baumannii ATCC 17978 plasmid pAB2, complete sequence' => 'dce',
'gi|126640115|ref|NC_009085.1| Acinetobacter baumannii ATCC 17978, complete genome' => 'dce',
'gi|169302972|ref|NC_010401.1| Acinetobacter baumannii AYE plasmid p1ABAYE, complete sequence' => 'dce',
'gi|169302980|ref|NC_010402.1| Acinetobacter baumannii AYE plasmid p2ABAYE, complete sequence' => 'dce',
'gi|169302992|ref|NC_010403.1| Acinetobacter baumannii AYE plasmid p4ABAYE, complete sequence' => 'dce',
'gi|169786889|ref|NC_010404.1| Acinetobacter baumannii AYE plasmid p3ABAYE, complete sequence' => 'dce',
'gi|169794206|ref|NC_010410.1| Acinetobacter baumannii AYE, complete genome' => 'dce',
'gi|169302963|ref|NC_010395.1| Acinetobacter baumannii SDF plasmid p1ABSDF, complete sequence' => 'dce',
'gi|169786833|ref|NC_010396.1| Acinetobacter baumannii SDF plasmid p2ABSDF, complete sequence' => 'dce',
'gi|169786864|ref|NC_010398.1| Acinetobacter baumannii SDF plasmid p3ABSDF, complete sequence' => 'dce',
'gi|169632029|ref|NC_010400.1| Acinetobacter baumannii SDF, complete genome' => 'dce',
'gi|299768250|ref|NC_014259.1| Acinetobacter sp. DR1 chromosome, complete genome' => 'dce',
'gi|82701135|ref|NC_007614.1| Nitrosospira multiformis ATCC 25196 chromosome, complete sequence' => 'dce',
'gi|82703893|ref|NC_007615.1| Nitrosospira multiformis ATCC 25196 plasmid 1, complete sequence' => 'dce',
'gi|82703911|ref|NC_007616.1| Nitrosospira multiformis ATCC 25196 plasmid 2, complete sequence' => 'dce',
'gi|82703928|ref|NC_007617.1| Nitrosospira multiformis ATCC 25196 plasmid 3, complete sequence' => 'dce',
'gi|124265193|ref|NC_008825.1| Methylibium petroleiphilum PM1 chromosome, complete genome' => 'bce',
'gi|124262546|ref|NC_008826.1| Methylibium petroleiphilum PM1 plasmid RPME01, complete sequence' => 'bce',
'gi|152979768|ref|NC_009659.1| Janthinobacterium sp. Marseille, complete genome' => 'dce',
'gi|134093294|ref|NC_009138.1| Herminiimonas arsenicoxydans chromosome, complete genome' => 'dce',
'gi|300309346|ref|NC_014323.1| Herbaspirillum seropedicae SmR1 chromosome, complete genome' => 'dce',
'gi|307069503|ref|NC_014497.1| Candidatus Zinderia insecticola CARI chromosome, complete genome' => 'dce',
'gi|116334902|ref|NC_008512.1| Candidatus Carsonella ruddii PV, complete genome' => 'bce',
'gi|305672698|ref|NC_014479.1| Bacillus subtilis subsp. spizizenii str. W23 chromosome, complete genome' => 'o',
'gi|49175990|ref|NC_000913.2| Escherichia coli str. K-12 substr. MG1655 chromosome, complete genome' => 'o',
'gi|190570478|ref|NC_010981.1| Wolbachia endosymbiont of Culex quinquefasciatus Pel, complete genome' => 'w',
'gi|42519920|ref|NC_002978.6| Wolbachia endosymbiont of Drosophila melanogaster, complete genome' => 'w',
'gi|58584261|ref|NC_006833.1| Wolbachia endosymbiont strain TRS of Brugia malayi, complete genome' => 'w',
'gi|225629872|ref|NC_012416.1| Wolbachia sp. wRi, complete genome' => 'w',);

#read in the blastn report and record qualifying reads
my(%blastnwreads,%blastndcereads,$in,$result,$hit,$hsp,$tothsplen);
$in = new Bio::SearchIO(-format => 'blast', -file   => $blastnrep);
$ctr=0; print STDERR "Reading blastn results..\n";
while($result = $in->next_result) {## $result is a Bio::Search::Result::ResultI compliant object
	if($result->no_hits_found()){next;}
	
	#get hit data
	if($result->num_hits>0){
		#debug
		if($debug){print STDERR "For ",$result->query_name,"\n";}
		@temp=();
		while($hit = $result->next_hit ) {# $hit is a Bio::Search::Hit::HitI compliant object
			$i=$hit->name().' '.$hit->description();#CHECK
			if($genomes{$i} eq 'w'){#if hit is a wolbachia species
				while($hsp = $hit->next_hsp()){# $hsp is Bio::Search::HSP::HSPI compliant object
					if($hsp->length('query') >= $mlen){#if length of read involved in HSP > match length
						#record name, ref genome strand/start/end/
						#recording both start and end since the inslen for each lib does not include read lengths
						push @temp,$hit->name(); push @temp,$hsp->strand('subject');
						push @temp,$hsp->start('subject'); push @temp,$hsp->end('subject');
					}
				}
			}
	    }
	    #debug
		if($debug){foreach $i (0..$#temp){if(($i%4==0)&&($i!=0)){print STDERR "\n"; } print STDERR "\t$temp[$i]";}}
	    #store in hash
		$blastnwreads{$result->query_name()}=[@temp];
	    $ctr++; if($ctr%100000 == 0){print STDERR "$ctr..";}
	    #debug
		if($debug){print STDERR "\n";}
	}
}
print STDERR "Read $ctr blastn results..\n";
if($debug){ $i=scalar (keys (%blastnwreads)); print STDERR "Length of blastnwreads hash: $i\n";}
#compute the qualifying reads
my(%qualwreads,%qualdcereads,%qualwPEcontamreads,$fwd,@fwdarr,$rev,@revarr);
$ctr=0; print STDERR "Searching for qualifying read pairs..";
if($debug){print STDERR "\n";}
while(($i,$j)=each(%blastnwreads)){
	if($i=~ /1$/){#forward read
		$fwd=$i; @fwdarr=@$j;
		$i=~ s/1$/2/;
		if(exists $blastnwreads{$i}){# get rev read data
			$rev=$i; $j=$blastnwreads{$i}; @revarr=@$j;
			#name,strand,start,end triplets in 1D array
			for($i=0;$i<@fwdarr;$i+=4){
				for($j=0;$j<@revarr;$j+=4){
					#if hit on same genome and diff strand
					if(($fwdarr[$i] eq $revarr[$j]) && ($fwdarr[$i+1] != $revarr[$j+1])){ 
						#fwd on pos and rev on comp strand for innies
						if (($fwdarr[$i+1] == 1) && ($revarr[$j+1] == -1) && ($fwdarr[$i+3]+$minins <= $revarr[$j+2]) 
							&& ($fwdarr[$i+3]+$maxins >= $revarr[$j+2])){
							$qualwreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
							$qualwreads{$rev}=$revarr[$j+2]-$fwdarr[$i+3];# rev read name and insert length
							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwreads{$rev}. fwd on pos and rev on comp strand for INNIES\n";}
							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT1;#done for this pair so get out
						}
#						#fwd on neg and rev on pos strand for innies
#						#decided to do this because lots of fwd.rev hits on comp/pos with proper insert size
#						elsif (($fwdarr[$i+1] == -1) && ($revarr[$j+1] == 1) && ($revarr[$j+3]+$minins <= $fwdarr[$i+2]) 
#							&& ($revarr[$j+3]+$maxins >= $fwdarr[$i+2])){
#							$qualwreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
#							$qualwreads{$rev}=$fwdarr[$i+2]-$revarr[$j+3];# rev read name and insert length
#							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwreads{$rev}. fwd on COMP and rev on POS strand for INNIES\n";}
#							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT1;#done for this pair so get out
#						}
						#fwd on comp and rev on pos for outies that were originally innie contaminants in a MP lib 
						#before REVCOMPING the MP lib 
						elsif (($mp) && ($fwdarr[$i+1] == -1) && ($revarr[$j+1] == 1) && ($fwdarr[$i+3]+$minins <= $revarr[$j+2])
							&& ($fwdarr[$i+3]+$maxins >= $revarr[$j+2])){
							$qualwPEcontamreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
							$qualwPEcontamreads{$rev}=$revarr[$j+2]-$fwdarr[$i+3];# rev read name and insert length
							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwPEcontamreads{$rev}. fwd on comp and rev on pos strand for OUTIES\n";}
							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT1;#done for this pair so get out
							
						}
#						#fwd on pos and rev on comp for outies that were originally innie/PE contaminants in a MP lib 
#						#before REVCOMPING the MP lib 
#						elsif (($mp) && ($fwdarr[$i+1] == 1) && ($revarr[$j+1] == -1) && ($fwdarr[$i+3]+$minins <= $revarr[$j+2])
#							&& ($fwdarr[$i+3]+$maxins >= $revarr[$j+2])){
#							$qualwPEcontamreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
#							$qualwPEcontamreads{$rev}=$revarr[$j+2]-$fwdarr[$i+3];# rev read name and insert length
#							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwPEcontamreads{$rev}. fwd on POS and rev on COMP strand for OUTIES\n";}
#							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT1;#done for this pair so get out
#							
#						}
						else{
							if($debug){
								print STDERR "\tRejected Fwd: $fwd and Rev: $rev due to insert length: ";
								print STDERR $revarr[$j+2]-$fwdarr[$i+3],"\tFwd strand: ",$fwdarr[$i+1],"\tRev strand: ",$revarr[$j+1],"\n";
#								if(($fwdarr[$i+1] == 1) && ($revarr[$j+1] == -1)){
#									print STDERR $revarr[$j+2]-$fwdarr[$i+3],"\tFwd strand: ",$fwdarr[$i+1],"\tRev strand: ",$revarr[$j+1],"\n";
#								}
#								elsif(($fwdarr[$i+1] == -1) && ($revarr[$j+1] == 1)){
#									print STDERR $fwdarr[$i+2]-$revarr[$j+3],"\tFwd strand: ",$fwdarr[$i+1],"\tRev strand: ",$revarr[$j+1],"\n";
#								}
							}
						}
					}
				}
			}
		}
		OUT1:
		#delete both reads from hash
		delete $blastnwreads{$rev};
		delete $blastnwreads{$fwd};
	}
	elsif($i=~ /2$/){#reverse read
		$rev=$i; @revarr=@$j;
		$i=~ s/2$/1/;
		if(exists $blastnwreads{$i}){# get fwd read data
			$fwd=$i; $j=$blastnwreads{$i}; @fwdarr=@$j;
			#name,strand,location triplets in array
			for($i=0;$i<@fwdarr;$i+=4){
				for($j=0;$j<@revarr;$j+=4){
					#if hit on same genome and diff strand
					if(($fwdarr[$i] eq $revarr[$j]) && ($fwdarr[$i+1] != $revarr[$j+1])){ 
						#fwd on pos and rev on comp strand for innies
						if (($fwdarr[$i+1] == 1) && ($revarr[$j+1] == -1) && ($fwdarr[$i+3]+$minins <= $revarr[$j+2]) 
							&& ($fwdarr[$i+3]+$maxins >= $revarr[$j+2])){
							$qualwreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
							$qualwreads{$rev}=$revarr[$j+2]-$fwdarr[$i+3];# rev read name and insert length
							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwreads{$rev}. fwd on pos and rev on comp strand for INNIES\n";}
							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT2;#done for this pair so get out
						}
#						#fwd on neg and rev on pos strand for innies
#						#decided to do this because lots of fwd.rev hits on comp/pos with proper insert size
#						elsif (($fwdarr[$i+1] == -1) && ($revarr[$j+1] == 1) && ($revarr[$j+3]+$minins <= $fwdarr[$i+2]) 
#							&& ($revarr[$j+3]+$maxins >= $fwdarr[$i+2])){
#							$qualwreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
#							$qualwreads{$rev}=$fwdarr[$i+2]-$revarr[$j+3];# rev read name and insert length
#							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwreads{$rev}. fwd on COMP and rev on POS strand for INNIES\n";}
#							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT2;#done for this pair so get out
#						}
						#fwd on comp and rev on pos for outies that were originally innie contaminants in a MP lib 
						#before REVCOMPING the MP lib 
						elsif (($mp) && ($fwdarr[$i+1] == -1) && ($revarr[$j+1] == 1) && ($fwdarr[$i+3]+$minins <= $revarr[$j+2])
							&& ($fwdarr[$i+3]+$maxins >= $revarr[$j+2])){
							$qualwPEcontamreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
							$qualwPEcontamreads{$rev}=$revarr[$j+2]-$fwdarr[$i+3];# rev read name and insert length
							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwPEcontamreads{$rev}. fwd on comp and rev on pos strand for OUTIES\n";}
							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT2;#done for this pair so get out
							
						}
#						#fwd on pos and rev on comp for outies that were originally innie/PE contaminants in a MP lib 
#						#before REVCOMPING the MP lib 
#						elsif (($mp) && ($fwdarr[$i+1] == 1) && ($revarr[$j+1] == -1) && ($fwdarr[$i+3]+$minins <= $revarr[$j+2])
#							&& ($fwdarr[$i+3]+$maxins >= $revarr[$j+2])){
#							$qualwPEcontamreads{$fwd}=$fwdarr[$i];#record fwd read name and genome it aligned to
#							$qualwPEcontamreads{$rev}=$revarr[$j+2]-$fwdarr[$i+3];# rev read name and insert length
#							if($debug){print STDERR "\tSelected Fwd: $fwd and Rev: $rev with insert length: $qualwPEcontamreads{$rev}. fwd on POS and rev on COMP strand for OUTIES\n";}
#							$ctr++; if($ctr%1000 == 0){print STDERR "$ctr..";} goto OUT1;#done for this pair so get out
#							
#						}
						else{
							if($debug){
								print STDERR "\tRejected Fwd: $fwd and Rev: $rev due to insert length: ";
								print STDERR $revarr[$j+2]-$fwdarr[$i+3],"\tFwd strand: ",$fwdarr[$i+1],"\tRev strand: ",$revarr[$j+1],"\n";
#								if(($fwdarr[$i+1] == 1) && ($revarr[$j+1] == -1)){
#									print STDERR $revarr[$j+2]-$fwdarr[$i+3],"\tFwd strand: ",$fwdarr[$i+1],"\tRev strand: ",$revarr[$j+1],"\n";
#								}
#								elsif(($fwdarr[$i+1] == -1) && ($revarr[$j+1] == 1)){
#									print STDERR $revarr[$j+2]-$fwdarr[$i+3],"\tFwd strand: ",$fwdarr[$i+1],"\tRev strand: ",$revarr[$j+1],"\n";
#								}
							}
						}
					}
				}
			}
		}
		OUT2:
		#delete both reads from hash
		delete $blastnwreads{$fwd};
		delete $blastnwreads{$rev};
	}
}
print STDERR "\nFound $ctr qualifying read pairs..\n";

#parse reads and write out the qualifying reads
my($fout,$rout,$seq,$uout,$fPEcontamout,$rPEcontamout);
$fout = Bio::SeqIO->new(-file=>">fwd.wol.${out}", -format=>$format);#NEED SEPARATE READS FILES??
$rout = Bio::SeqIO->new(-file=>">rev.wol.${out}", -format=>$format);
$uout = Bio::SeqIO->new(-file=>">unmapped.${out}", -format=>$format);

if($mp){
	$fPEcontamout = Bio::SeqIO->new(-file=>">fwd.PEcontam.wol.${out}", -format=>$format);#NEED SEPARATE READS FILES??
	$rPEcontamout = Bio::SeqIO->new(-file=>">rev.PEcontam.wol.${out}", -format=>$format);
}
print XLS "\n\nPair No\tFwd read\tRev read\tGenome\tInsert size\n";
$i=$ctr=$l=0; print STDERR "Writing qualifying and unmapped reads..";
$in = Bio::SeqIO->new(-file=>"<$pairedreads", -format=>$format);
while ($seq = $in->next_seq()){
	if(exists $qualwreads{$seq->display_id()}){
		if($seq->display_id()=~ /1$/){
			if($format eq 'fastq'){ $fout->write_fastq($seq);}
			else{$fout->write_seq($seq);}
			$j=$qualwreads{$seq->display_id()};
		}
		elsif($seq->display_id()=~ /2$/){
			if($format eq 'fastq'){ $rout->write_fastq($seq);}
			else{$rout->write_seq($seq);}
			$k=$seq->display_id(); $k=~ s/2$/1/;
			#writing on rev since genome recorded from fwd
			print XLS (($ctr+1)/2),"\t",$k,"\t",$seq->display_id(),"\t",$j,"\t",$qualwreads{$seq->display_id()},"\n";
		}
		$ctr++; if($ctr%10000 == 0){print STDERR "$ctr..";}
	}
	elsif(($mp) && (exists $qualwPEcontamreads{$seq->display_id()})){
		if($seq->display_id()=~ /1$/){
			if($format eq 'fastq'){ $fPEcontamout->write_fastq($seq);}
			else{$fPEcontamout->write_seq($seq);}
			$j=$qualwPEcontamreads{$seq->display_id()};
		}
		elsif($seq->display_id()=~ /2$/){
			if($format eq 'fastq'){ $rPEcontamout->write_fastq($seq);}
			else{$rPEcontamout->write_seq($seq);}
			$k=$seq->display_id(); $k=~ s/2$/1/;
			#writing on rev since genome recorded from fwd
			print XLS (($ctr+1)/2)," PEcontam\t",$k,"\t",$seq->display_id(),"\t",$j,"\t",$qualwPEcontamreads{$seq->display_id()},"\n";
		}
		$ctr++; $l++; if($ctr%10000 == 0){print STDERR "$ctr..";}
	}
	else{
		if($format eq 'fastq'){ $uout->write_fastq($seq);}
		else{$uout->write_seq($seq);}
		$i++;
	}
}
close(XLS);
if($mp){print STDERR "\nWrote $ctr qualifying MP reads, $l qualifying innie MP reads and $i unmapped reads\n";}
else{print STDERR "\nWrote $ctr qualifying PE reads and $i unmapped reads\n";}

my($user_t,$system_t,$cuser_t,$csystem_t); ($user_t,$system_t,$cuser_t,$csystem_t) = times;
print STDERR "System time for process: ",sprintf("%.3f",$system_t/3600)," hrs\n"; print STDERR "User time for process: ",sprintf("%.3f",$user_t/3600)," hrs\n";

exit;