#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:r:s:l:u:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{r}=1 unless defined($opts{r});
$opts{s}=18 unless defined($opts{s});
$opts{l}=26 unless defined($opts{l});
$opts{u}="rpm" unless defined($opts{u});

###############################################################################
#Main text.
###############################################################################
my @filelist=glob "$opts{i}";
my %database;
# Befere filter abundance
$t=0; # record total reads
$u=0; # record unique count
foreach (@filelist){
	open INPUT,"<$_";
	$basename=basename($_);
#	open OUTPUT,">$opts{o}$basename";
	print $_,"\n";
	while (<INPUT>){
		if ($_=~/^([AGCT]{$opts{s},$opts{l}})\t(\d+)/){
				$database{$1}+=$2;
				$t+=$2;
#			if ($2>=$opts{r}) {
#				print OUTPUT "$1\t$2\n";
#			}
		}
	}
	close INPUT;
#	close OUTPUT;
}
if ($opts{u} eq "rpm") {
	$opts{rpm}=int($opts{r}*$t/1000000);
}else{
	$opts{rpm}=$opts{r};
}

open OUTPUT,">$opts{o}";
# After filter abundance
$t1=0; # record total reads
$u1=0; # record unique count
foreach (keys %database){
	$u++;
	if ($database{$_} > $opts{rpm}) {
		$t1+=$database{$_};
		$u1++;
		print OUTPUT "$_\t$database{$_}\n" ;
	}
}
close OUTPUT;
print "Threshold of length is $opts{s} to $opts{l}, abundance is $opts{r} RPM and $opts{rpm} reads in average.\n";
print "Total reads $t, unique reads $u\nFilter reads $t1, unique reads $u1\n";
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
        qq!
Usage:    merge2.pl -i inpute_file_list -o output_dir
Function: Merge all reads count into reads count, filter 18-26 length and abundance more than 1
Command:  -i inpute file list (Must)
          -o output directory
          -r reads filter threshold, default=1, select>1 reads/rpm, length 18-26
		  -u unit, filter reads or rpm, default rpm
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/7/11
Notes:    
\n!
    )
}
