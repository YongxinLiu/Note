#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:r:', \%opts );
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
	open OUTPUT,">$opts{o}$basename";
	print $_,"\n";
	while (<INPUT>){
		if ($_=~/^([AGCT]{20,24})\t(\d+)/){
			if ($2>=$opts{r}) {
				$database{$1}+=$2;
				print OUTPUT "$1\t$2\n";
			}
		}
	}
	close INPUT;
	close OUTPUT;
}
open OUTPUT,">$opts{o}merge";
foreach (keys %database){
	print OUTPUT "$_\t$database{$_}\n";
}
close OUTPUT;
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
Usage:    merge2.pl -i inpute_file_list -o output_dir
Function:
Command:  -i inpute file list (Must)
          -o output directory
          -r reads filter threshold, default=2, select>=2 reads, length 20-24
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-09-03
Notes:    
\n/
    )
}
