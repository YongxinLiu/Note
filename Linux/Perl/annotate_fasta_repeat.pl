#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my %database; #database in hash
#SW	perc	perc	perc	query	position	in	query	matching	repeat	position	in	repeat
#score	div.	del.	ins.	sequence	begin	end	(left)	repeat	class/family	begin	end	(left)	ID
#
#912	0.0	1.0	0.0	0_2093	1	101	(0)	C	ANGELA6_TM-LTR	LTR/Copia	(981)	740	639	1
#806	5.0	0.0	0.0	10_736	1	101	(0)	+	Gypsy-12_TA-LTR	LTR/Gypsy	1384	1484	(37)	2
#814	3.0	1.0	1.0	11_725	1	101	(0)	C	ANGELA6_TM-LTR	LTR/Copia	(967)	754	654	3
#425	10.4	3.0	5.0	12_705	1	101	(0)	C	WIS21A_TA-int	LTR/Copia	(4243)	871	773	4
#737	10.9	0.0	0.0	16_655	1	101	(0)	C	SSU-rRNA_Ath	rRNA	(365)	1537	1437	5	
#369	17.5	1.2	4.9	17_655	1	84	(17)	C	WIS21A_TA-int	LTR/Copia	(4355)	759	679	6
#874	3.0	0.0	0.0	19_629	1	101	(0)	C	EnSpm-4_HV	DNA/CMC-EnSpm	(9196)	2539	2439	7
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[4]}="$tmp[10]\|$tmp[9]";
}

#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
close DATABASE;

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
#>0_2093
#CCTACATATTCTACGAAGATCTTTATCGGTCAGACCGCATAACAACATACGTTGTTCCCTTTGTCATCGGTATGTTACTTGCCCGAGATTCGATCGTCGGT
#>1_1782
#ATATAAGTGATGTGATATGGTCAAGACATGATGCTAAATTTTATTGTATGAGATGATCATGTTTTGTAACCGAGTTATCGGCAACTGGCAGGAGCCATATG
#>2_1084
#CTTTGTTATTTGGTTGCAGGGTTGCTTGAGAGAGACCATCTTCATCCTACGCCTCCTACGGATTGATAAACCTTAGGTCATCCACTTGAGGGAAATTTGCT
#>3_1070
#GACCTTCTTCATGGACACCACGTGCATCATCTACTTCCACACGTTCTGCTTCACCGGCGCCTCCGTCGGCTGCCACCTCCAACCGTGATACAATAAAGCAG
#>4_991
#CAACTAAACATTCATACCTGTCAGCATCTTCAGGAGCCATGTATTGCGTCTCATGATCAGAATCGTCTCCACCGTGTTCGTCACGTGTAATAAGAGCCAAA
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	$_=~/>([\w\_]+)/;
	$seq=<INPUT>;
	chomp($seq);
	print OUTPUT "$1\t$seq\t$database{$1}\n";
	
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
Usage:    annotate_fasta_repeat.pl -i fasta -o table -d repeat output (sep  by table) -h header num
Function: annotate fasta sequence by repeat
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}