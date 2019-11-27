#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:c:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{c}=0 unless defined($opts{c});

my (@tmp1);
open DATABASE,"<$opts{d}";
while (<DATABASE>) {
	chomp;
	my @tmp=split/\s+/;
	print $tmp[0],"\n";
	push @tmp1,$tmp[$opts{c}];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
my $i=0;
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	foreach $x(@tmp1) {
		if ($_=~m/$x/g) {
			print OUTPUT $_;
			last;
		}
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
Usage:    search_line.pl -i line_file -d sequence or gene ID -o line_file
Function: search sequence or gene ID from a dataset
Command:  -i line_file (Must)
          -d GI (Must)
          -o searching result
          -c column, default=0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-09-15
Notes:    Add column select GI column in database
\n/
    )
}
