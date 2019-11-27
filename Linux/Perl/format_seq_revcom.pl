#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i});
my $start_time=time;
$opts{o}="$opts{i}.out" unless defined($opts{o});
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";

open OUTPUT,">$opts{o}";
my %seq;
while (<INPUT>) {
	chomp;
	$seq=$_;
	$seq=~tr/ACGT/TGCA/;
	$seq=reverse($seq);
	print "$seq\n";
	print OUTPUT "$seq\n";
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
        qq!
Usage:    format_seq_revcom.pl -i inpute_file -o output_file
Function: fomrat sequence into reverse complement by each line
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/12/12
Notes:    
\n!
    )
}