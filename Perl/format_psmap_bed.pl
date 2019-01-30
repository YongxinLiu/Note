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
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#id&reads0	chr1	strand2		start3	end4 	sequences5
#21719662&1      1       +       656     677     TGGATATAAACTATTTTTGGCT
open OUTPUT,">$opts{o}";
#chrom Start End name score strand; socre was replace by count
#chr22 1000 5000 cloneA 960 + 
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	my @tmp1=split(/&/,$tmp[0]);
	print OUTPUT "$tmp[1]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp1[1]\t$tmp[2]\n";
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
Usage:    format_psmap_bed.pl -i inpute_file -o output_file
Function: only transformat, not need database for loci filter
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-23
Notes:    Only transform format
\n/
    )
}
