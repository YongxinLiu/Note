#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:r:m:n:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{r}=2 unless defined($opts{r});
$opts{m}=18 unless defined($opts{m});
$opts{n}=26 unless defined($opts{n});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>){
	if ($_=~/^([AGCT]{$opts{m},$opts{n}})\t(\d+)/){
		if ($2>=$opts{r}) {
			print OUTPUT "$1\t$2\n";
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
Usage:    filter_sRNA.pl -i inpute_file -o output_file
Function:
Command:  -i inpute file name
          -o output file name
          -r reads filter threshold, default=2, select>=2 reads
          -m min reads length, default=18
          -n max reads length, default=26
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-02
Notes:    
\n/
    )
}
