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
$opts{t}=1 unless defined($opts{t});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my ($total1,$total2)=(0,0);
my (%sample1,%sample2,%all);
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$sample2{$tmp[0]}=$tmp[1];
	$all{$tmp[0]}+=$tmp[1];
	$total2+=$tmp[1];
}
close DATABASE;

open DATABASE,"<$opts{i}";
#database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$sample1{$tmp[0]}=$tmp[1];
	$all{$tmp[0]}+=$tmp[1];
	$total1+=$tmp[1];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
foreach (keys %all) {
	next if $all{$_}<$opts{t};
	$sample1{$_}=0 unless defined($sample1{$_});
	$sample2{$_}=0 unless defined($sample2{$_});
	$sample1{$_}=$sample1{$_}/$total1*1000000;#read to RPM
	$sample2{$_}=$sample2{$_}/$total2*1000000;
	next if ($sample1{$_}<$opts{t} && $sample2{$_}<$opts{t});
	printf OUTPUT "$_\t%.2f\t%.2f\n",$sample1{$_},$sample2{$_};
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
Usage:    compare_replicate_sRNA1.pl -i sample1 -d sample2 -o output_file -t threshold of reads
Function: compare two samples replicate relation, filter RPM>threshold, and output RPM
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -t threshold RPM, default >= 1
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-05-03
Notes:    Change reads to RPM as threshold
\n/
    )
}
