#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
$opts{d}=0.3 unless defined($opts{d});
print "Select flank scale $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

##database in array
#my (@tmp1,@tmp2);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;


#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr0	1				2		start3	end4	RPM5	strand6	7		ID8								left9	
#Chr1    omicslab        SRC     1       201     9.06    +       .       ID=ath-SRC1;length_type=24-nt   ath-rep1|Chr1|1|107|-|107|Repbase13|DNA|ATREP18;
open OUTPUT,">$opts{o}";
#ID0				overlap1(overlap+ID+description)			
#ath-SRC257      317|AT1G05780|Vacuolar ATPase assembly integral membrane protein VMA21-like domain;1|ath-rep747|A-rich 
while(<INPUT>){
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=(\w+\-\w+)/; #match ID
	print OUTPUT "$1\t";
	if ($tmp[9]=~/intergenic/) {#output intergenic
		print OUTPUT "intergenic\n";
		next;
	}
	my ($type,$gi,$description,%annotation);
	$overlap=0;
	my @tmp1=split/\;/,$tmp[9]; #split all annotation
	$pointer=5;
	foreach $anno(@tmp1) {
#		print $anno,"\n";
		my @tmp2=split/\|/,$anno; #split one annotation detail
#		print $tmp2[0],"\t",$tmp2[8],"\n";
		$tmp2[8]=$tmp2[7] unless defined($tmp2[8]);
		$annotation{$tmp2[0]}=$tmp2[8];
		next unless $tmp2[5]>$overlap;
		if ($anno=~/RNA$|RNA\||\|tRNA/ && $tmp2[5]>$opts{d}*($tmp[4]-$tmp[3]+1) && $tmp2[5]>20) { #match all the TAIR annotated RNA(RNA$) and miRNA, overlap>30%,
			$type="kncRNA"; #known non-coding RNA
			$overlap=$tmp2[5];
			$gi=$tmp2[0];
			$description=$tmp2[8];
			$pointer=1;
		}elsif ($anno=~/Repbase|mtec|transposable|pseudogene/ && $tmp2[5]>$opts{d}*($tmp[4]-$tmp[3]+1) && $tmp2[5]>20) { #repeat, and overlap>30%
			next if $pointer<2;
			$type="repeat";
			$overlap=$tmp2[5];
			$gi=$tmp2[0];
			$description=$tmp2[8];
			$pointer=2;
		}elsif ($anno=~/intron/ && $tmp2[5]>$opts{d}*($tmp[4]-$tmp[3]+1) && $tmp2[5]>20) { #intron, and overlap>30%
			next if $pointer<3;
			$type="intron";
			$overlap=$tmp2[5];
			$description=$tmp2[8];
#			print $description,"\n";
			$gi=$tmp2[0];
			$pointer=3;
		}elsif ($anno=~/TAIR10|ensembl/ && $tmp2[5]>$opts{d}*($tmp[4]-$tmp[3]+1) && $tmp2[5]>20) { #gene, and overlap>30%
			next if $pointer<4;
			$type="exon";
			$overlap=$tmp2[5];
			$gi=$tmp2[0];
			$description=$tmp2[8];
			$pointer=4;
		}else{
			$pointer=5;
		}
	}
	if ($pointer==3) {
		$description=$annotation{$description};
	}
	if ($pointer==5) {
		print OUTPUT "intergenic\n";
	}else{
		print OUTPUT "$type\|$overlap\|$gi\|$description\n";
	}
}

close INPUT;
close OUTPUT;

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq/
Usage:    Annotate_SRCs_overlap_simplify1.1.pl -i inpute_file -o output_file
Function: simplify overlap annotate region, add scale select
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d scale overlap percentage, default is more than 30%
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-15
Notes:    Imrpove gene region annotation, all the transposon and pseudogene as repeat
\n/
    )
}