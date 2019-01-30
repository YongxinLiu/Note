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
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#ath-SRC330	1.05	6.26	6.83
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$ave=($tmp[1]+$tmp[2]+$tmp[3])/3;#calculate average
	if ($tmp[1]>=0.5*$ave && $tmp[1]<=2*$ave && $tmp[2]>=0.5*$ave && $tmp[2]<=2*$ave && $tmp[3]>=0.5*$ave && $tmp[3]<=2*$ave) {
		print OUTPUT "$_\t0\n";
	}elsif ($tmp[1]>=$tmp[2] && $tmp[1]>=$tmp[3]) {
		print OUTPUT "$_\t1\n";
	}elsif ($tmp[2]>=$tmp[1] && $tmp[2]>=$tmp[3]){
		print OUTPUT "$_\t2\n";
	}elsif ($tmp[3]>=$tmp[2] && $tmp[3]>=$tmp[1]){
		print OUTPUT "$_\t3\n";
	}
	
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
Usage:    group_SRC.pl -i SRCs_RPM -o output_file
Function: group SRCs RPM in different groups, 0 is not significant changed, 1 is descreased, 2 is middle highest, 3 is increased.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-12
Notes:    
\n/
    )
}