#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:r:', \%opts );
&usage unless ( exists $opts{i} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
$opts{r}=0 unless defined($opts{r});

###############################################################################
#Read the database in memory(opt)
###############################################################################
# open a list file
my %list;
my @filelist=glob "$opts{i}";
foreach $file(@filelist){
	open DATABASE,"<$file";
	$file=basename($file);
	$file=~/^(.+)\.seq/;
#	$file=~/^([\w\-]+)(.+)\.seq/;
	$file=$1;
	while (<DATABASE>) {
		chomp;
		$list{$file}.=$_;
	}
	close DATABASE;
}

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
foreach  (sort keys %list) {
	print OUTPUT ">$_\n";
	$seq=$list{$_};
	if ($opts{r}==0) {
		print OUTPUT "$seq\n";
	}else{
		$seq=reverse($seq);
		$seq=~tr/ACGT/TGCA/;
		print OUTPUT "$seq\n";
	}
}
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
Usage:    format_seq2fasta.pl -i filelist -o output_file.fasta
Function: format batch sanger sequencing into fasta
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -r revcom, default is 0 off, 1 is on
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-08-04
Notes:    
\n/
    )
}
