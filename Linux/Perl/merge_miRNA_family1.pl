#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

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

#open a list file
#my @filelist=glob "$opts{i}";
#foreach (@filelist){
#	print "$_\n";
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#zma-miR171j-3p  5.19    11.22   7.81    5.69    4.02    1.46    5.62    3.27    87.73
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	$tmp=<INPUT>;
	$opts{h}--;
	print OUTPUT $tmp;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$length=$#tmp unless defined($length);
	$tmp[0]=~/([\w-]{7}\d+)(\w)?(\-\dp)?/;
	$id=$1.$3;
#	print "$1\t$2\t$3\n";
	for (1..$length) {
		$family{$id}[$_]+=$tmp[$_];
	}
}

foreach $id(sort keys %family) {
	print OUTPUT "$id";
	for (1..$length) {
		print OUTPUT "\t$family{$id}[$_]";
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
Usage:    merge_miRNA_family1.pl -i miRNA_normReads -o miRNA_family_reads
Function: merge each miRNA member reads\/RPM to family, split 3p and 5p
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-21
Notes:    
\n/
    )
}