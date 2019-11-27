#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:s:e:b:c:l:', \%opts );
&usage unless ( exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{b}=200 unless defined($opts{b});
$opts{l}=100 unless defined($opts{l});



###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
my $count=0;
#9	55456000	56181000	sDic15_725kb	725
while ($opts{s}+$opts{l}<$opts{e}) {
	$count++;
	$end=$opts{s}+$opts{b};
	$start=$opts{s};
	print OUTPUT "$opts{c}\t$start\t$end\t$opts{c}\_$start\_$end\t\.\t\+\n";
	$opts{s}+=$opts{l};
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
Usage:    write_bed_bin.pl -c chr -s start -e end -o output_file -b bin_size
Function: write a constitute bed file from start to end, step by bin
Command:  -s start (Must)
          -e end (Must)
          -c chr (Must)
          -b bin_size, default=200
          -l slid step, default=100
          -o output file name (Must)
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-08
Notes:    bed文件的起始要比真实位置小1，终点一致。如基因位置3:11-20，写成bed应3:10-20
\n/
    )
}