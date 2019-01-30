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
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
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
#>AGGGAGAGGAAHTTCGGG	1-AGGGAGAGGAAHTTCGGG,BestGuess:Trl/dmmpmm(Bigfoot)/fly(0.638)	13.955561	-105.500604	0	T:22.0(0.63%),B:1.7(0.00%),P:1e-45
#0.672	0.128	0.199	0.001
#0.085	0.001	0.829	0.085
open OUTPUT,">$opts{o}";
#MEME version 4
#
#ALPHABET= ACGT
#
#strands: + -
#
#Background letter frequencies
#A 0.264 C 0.233 G 0.233 T 0.264
#
#MOTIF crp
#letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009 
# 0.000000  0.176471  0.000000  0.823529 
# 0.000000  0.058824  0.647059  0.294118 
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my ($id,$number,$pvalue);
while (<INPUT>) {
	chomp;
	if (/>/) {
		my @tmp=split/\s+/;
#		print $tmp[5],"\n";
		$tmp[5]=~/T:(\d+).*P:([\w\-]+)/;
#		print "$1\t$2\n";
		$number=$1;
		$pvalue=$2;
		$tmp[0]=~/>(\w+)/;
		$id=$1;
	}elsif (/\d/){
		push @matrix,$_;
	}
}
#		print $#matrix,"\n";
$len=$#matrix+1;
print OUTPUT "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.264 C 0.233 G 0.233 T 0.264\n\nMOTIF $id\nletter-probability matrix: alength= 4 w= $len nsites= $number E= $pvalue\n";
foreach  (@matrix) {
	@tmp=split/\s+/;
	print OUTPUT " $tmp[0]  $tmp[1]  $tmp[2]  $tmp[3] \n";
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
Usage:    format_homer2meme.pl -i homer_motif1.motif -o meme_motif.txt
Function: Format homer motif matrix to meme format, and used by fimo
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