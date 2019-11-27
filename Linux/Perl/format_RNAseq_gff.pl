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
#test_id0 locus1   FPKM_NL2 FPKM_SL3 FPKM_NR4 FPKM_SR5 P_NL_SL6 P_NR_SR7 P_NL_NR8 P_SL_SR9 annotation10      biological_process11
#XLOC_000001     Chr1:3630-5899  5.00205 5.27591 23.5308 38.8806 0.0336365696776894      6.03355050428681e-19    2.59784606653356e-123   1.29806438263668e-208   NAC domain containing protein 1;        multicellular organismal development, DNA-dependent, transport, developmental processes, transcription, regulation of transcription; 
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
open OUTPUT,">$opts{o}";
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#9	ensembl	chromosome	1	156750706	.	.	.	ID=9;Name=chromosome:AGPv2:9:1:156750706:1
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[1]=~/(\w+)\:(\d+)\-(\d+)/;
	print OUTPUT "$1\txiuying\tmRNA\t$2\t$3\t$tmp[2]\|$tmp[3]\|$tmp[4]\|$tmp[5]\|$tmp[6]\|$tmp[7]\|$tmp[8]\|$tmp[9]\t+\t\.\tID=$tmp[0];Anno=$tmp[10]GO=$tmp[11]\n";
#	print OUTPUT "$1\txiuying\t$2\t$3\t\.\t+\.\tID=$tmp[0]\n";
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
Usage:    format_RNAseq_gff.pl -i inpute_file -o output_file
Function: Format xiuying RNA-seq result table to gff
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -h header, default=0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-22
Notes:    
\n/
    )
}