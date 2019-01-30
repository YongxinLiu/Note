#!/usr/bin/perl
#use strict;
#use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Bio::SeqIO;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
open DATABASE,"<$opts{d}";
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split('\t',$_);
#	print $tmp[0],"\n",$tmp[6],"\n",$tmp[8],"\n";
#	$database{$tmp[0]}="$tmp[1]\|$tmp[2]\|$tmp[3]\|$tmp[4]";
	$database{$tmp[0]}="$tmp[2]";
}
close DATABASE;

while (<INPUT>) {
	my @tmp=split('\t',$_);
	$_=~s/Note\=[^\t]*\t/$database{$tmp[4]}\t/;
	print OUTPUT $_;
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
Function:
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-09-20
Notes:    
\n/
    )
}