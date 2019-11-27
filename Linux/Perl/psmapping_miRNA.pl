#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

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

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
my %rpkm;
#0			1
#76559&2 nta-MIR6020b MI0020235 Nicotiana tabacum miR6020b stem-loop     +       21      41      AAATGTTCTTCGAGTATCTTC
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	my @tmp1=split/&/,$tmp[0];
	my @tmp2=split/\s+/,$tmp[1];
	$rpkm{$tmp2[0]}+=$tmp1[1];
}
foreach  (sort keys %rpkm) {
	$rpm=$rpkm{$_}/$opts{t}*1000000;
	printf OUTPUT "$_\t$rpkm{$_}\t%.2f\n",$rpm;
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
Usage:    template.pl -i inpute_file -o output_file
Function: get rpm from psRobot_map
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -t total reads
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-05-04
Notes:    
\n/
    )
}