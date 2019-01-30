#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:l:r:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{l}=26 unless defined($opts{l});
$opts{r}=0 unless defined($opts{r});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#@SRR488770.1 GSM918104_r1:5:1:564:495 length=36
#GATGACTGGGCACACATGCTGGCATCGTATGCCGTC
#+SRR488770.1 GSM918104_r1:5:1:564:495 length=36
#IIIIIIIIIIIIIIIIIIIIIIIIII?I3I(9I/=1
open OUTPUT,">$opts{o}";
my %reads;
my $i=0;
while (<INPUT>) {
	chomp;
	if (/(^[AGCT]+$)/) {
		next if length($1)>$opts{l};
		$reads{$1}++;
	}
}
my $j=0;
foreach (keys %reads) {
	next if $reads{$_}<=$opts{r};
	$j++;
	$i+=$reads{$_};
	print OUTPUT "$_\t$reads{$_}\n";
}

print "Total reads $i\nUnique reads $j\n";

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
Usage:    format_fastq_sRNA.pl -i fastq -o sRNA.txt -d database	-h header num
Function: format fastq to reads and count, filter reads content N
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
          -l max length, default 26
          -r reads min, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-09-19
Notes:    
\n/
    )
}
