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
open DATABASE,"<$opts{d}";
my %database; #database in hash
while (<DATABASE>) {
	#LTR/Gypsy       TE
	#rRNA    TE
	#Satellite       SSR
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[1];
}

close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#0_25498 CATGAAACAAGTTTAAAGGACAAGCTGTCCAAAACAGCAACTAGACAACTTGCAAGATGACAACAACATCACACTACCAAGTTAAAGGTCTCTCATCATTA   
#1_9033  AGGAGAATGATGTGATTGACTTTACCCATTCCGTTAGACTAGCACTTGATCATTTTAGTTTACTGCTACTGCTTTCTTCATGACTTATACATGATCCTATG   LTR/Copia|ANGELA6_TM-LTR
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (!defined($tmp[2])) {
		print OUTPUT "$tmp[0]\t$tmp[1]\tUnknown\t\n";
		next;
	}
	my @tmp1=split(/\|/,$tmp[2]);
	if (defined($tmp1[0])) {
		print OUTPUT "$tmp[0]\t$tmp[1]\t$database{$tmp1[0]}\t$tmp[2]\n";
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
Usage:    classify_repeatmasker.pl -i AADD_us_top6000.repeat -o AADD_us_top6000.repeat.class -d repeat_group.txt -h header num
Function: Classify repeatmasker into three group, mainly TE, SSR and Unknown
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-11-04
Notes:    
\n/
    )
}