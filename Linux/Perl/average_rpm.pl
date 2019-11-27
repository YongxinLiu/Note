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
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
##############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	for  (1..$#tmp) {
		$database{$tmp[0]}[$_]=$tmp[$_];
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#ID				RPM		reads	copies
#zma-SRC1        2.19    19      1.31
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT $tmp[0];
	for  (1..$#tmp) {
		$average=($tmp[$_]+$database{$tmp[0]}[$_])/2;
		print OUTPUT "\t$average";
	}
	print OUTPUT "\n";
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
Usage:    average_rpm.pl -i a.rpm -d b.rpm -o average.rpm
Function: get the average of two samples
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-22
Notes:    
\n/
    )
}