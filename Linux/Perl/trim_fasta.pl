#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:s:l:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{s}=0 unless defined($opts{s});
$opts{l}=1000 unless defined($opts{l});

###############################################################################
#Read the database in memory(opt)
###############################################################################

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";

while (<INPUT>) {
	chomp;
	if (/>/){
		print OUTPUT "$_\n";
	}else{
		$tmp=substr($_,$opts{s},$opts{l});
		print OUTPUT "$tmp\n";
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
Usage:    trim_fasta.pl -i gene.fa -d 1.txt -o gene1.fa
Function: trim fasta to same length
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -s default 0
          -l length, default 1000
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-07-29
Notes:    
\n/
    )
}