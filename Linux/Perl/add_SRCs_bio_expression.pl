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
my %list;
my @filelist=glob "$opts{d}";
foreach $file(@filelist){
	$filename=basename($file);
	next unless $filename=~/([A-Za-z]+\d+)_([A-Za-z]+)\.diff/;
	print "$1\t$2\t$filename\n";
	$mutant="$1\-$2";
	open DATABASE,"<$file";
	while (<DATABASE>) {
#ID				RPM1	RPM2	Fold3	P-value4			FDR5
#ath-SRC10444    6.22    0.66    4.35    0.000000032383  0.00000029
	my @tmp=split/\t/;
		$list{$tmp[0]}{$mutant}="$tmp[1]\|$tmp[2]\|$tmp[3]";
	}
	close DATABASE;
}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT "$tmp[0]";
	foreach  $id(1..$#tmp) {
		if (defined($list{$tmp[0]}{$tmp[$id]})) {
			print OUTPUT "\t$tmp[$id]\|$list{$tmp[0]}{$tmp[$id]}";
		}else{
			print OUTPUT "\t$tmp[$id]";
		}
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
Usage:    add_SRCs_bio_expression.p -i inpute_file -d expression_list  -o output_file
Function: Add all the SRCs biogenesis expression
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-06-07
Notes:    
\n/
    )
}