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
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#@HWI-ST833:261:3:1101:1359:2128#0/1
#GACGGCTTCTCCTAGGCGCGGGGAGGGAGGCGGCCTCCGGTCTGGTGGAGCATGGCATCGATCTCGTCGCCATCTTCCATCTCAAGCTCATCAGGGGTCTG
#+
#@<@FFFFFDFFHBHIIDFIEHGG/@BBD<<?>?150ABBB&)7C<+<?+2?A9@BB9:A7>AABC?B@BBB???C@CCCDACCA@>CCC@@>C:>2<9599
	print OUTPUT $_;
	$seq=<INPUT>;
	chomp($seq);
	$revcom=reverse($seq);
	$revcom=~tr/AGCT/TCGA/;
	print OUTPUT "$revcom\n";
	$blank=<INPUT>;
	print OUTPUT $blank;
	$value=<INPUT>;
	chomp($value);
	$rev=reverse($value);
	print OUTPUT "$rev\n";
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
Usage:    format_fastq_revcom.pl -i fastq -o fastq
Function: reverse complement
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}