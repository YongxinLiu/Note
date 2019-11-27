#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:a:p:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################
my @filelist=glob "$opts{i}";
`mkdir fastq`;
`mkdir fasta`;
`mkdir clipper`;
for (0..$#filelist){
	my $basename=basename($filelist[$_]);
	my $output=$opts{o}.$basename;
	print "$filelist[$_]\tProgress in %",$_/@filelist*100,"!!!\n";
	print "fastq-dump $opts{i}$basename -Z >fastq/$basename\n";
	`fastq-dump $opts{i}$basename -Z >fastq/$basename`;
	print "format_fastq_fasta.pl -i fastq/$basename -o fasta/$basename\n";
	`format_fastq_fasta.pl -i fastq/$basename -o fasta/$basename`;
	#fastx_clipper -a TCGTATGCCGTC -l 15 -c -i GSM562942_9311s.fasta -o GSM562942_9311s.fasta.clipper #remove adaptor
	#print "$opts{p} -a $opts{a} -l 15 -c -i $filelist[$_] -o $output\n";
	print "fastx_clipper -a $opts{a} -l 18 -c -i fasta/$basename  -o clipper/$basename\n";
	`fastx_clipper -a $opts{a} -l 18 -c -i fasta/$basename  -o clipper/$basename`;
	print "format_fasta_table.pl -i clipper/$basename -o $opts{o}$basename\n";
	`format_fasta_table.pl -i clipper/$basename -o $opts{o}$basename`;
}
print scalar @filelist," files have been treated.\n";


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
Usage:    batch_sra_sRNA.pl -i sra_dir -o sRNA -a adaptor
Function: format SRA to sRNA format, based on fastx clipper adaptor
Command:  -i inpute_file_list (Must)
          -o output_file_directory
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-10-12
Notes:    
\n/
    )
}
