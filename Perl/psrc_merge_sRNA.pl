#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{r}=2 unless defined($opts{r});

###############################################################################
#Main text.
###############################################################################
my @filelist=glob "$opts{i}";
my %database;
foreach (@filelist){
	open INPUT,"<$_";
	$basename=basename($_);
	print $_,"\n";
	while (<INPUT>){
		if ($_=~/^([AGCT]+)\t(\d+)/){
				$database{$1}+=$2;
		}
	}
	close INPUT;
}
open OUTPUT,">$opts{o}";
foreach (keys %database){
	print OUTPUT "$_\t$database{$_}\n";
}
close OUTPUT;
print scalar @filelist," files have been merged.\n";


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
Usage:    psrc_merge_sRNA.pl -i inpute_file_list -o output_dir
Function: merge all the sample into one unique reads file
Command:  -i inpute file list (Must)
          -o output directory
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-02
Notes:    
\n/
    )
}
