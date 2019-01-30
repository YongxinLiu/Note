#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{h}=0 unless defined($opts{h});


##############################################################################
#Main text
###############################################################################
$t=$opts{d}/$opts{t};
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	if (rand()<$t) {
		print OUTPUT $_;
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
Usage:    sampling_line.pl.pl -i linefile -o sampled_file -d number -t total
Function: sampling about D line of total T line file
Command:  -i linefile
          -o output file name
          -d output line number
          -t total line
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-10-09
Notes:    
\n/
    )
}