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
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{c}=0 unless defined($opts{c});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my %database;
while (<DATABASE>) {
	chomp;
	if (/>/) {
		$_=~/>([^\s]+)/;
		$database{$1}=1;
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#>1
#TCGGACCAGGCTTCATTCCCC
open OUTPUT,">$opts{o}";
$pointer=0;
my %seq;
while (<INPUT>) {
	if (/>/){
		$_=~/>([^\s]+)/;
		if (defined($database{$1})){
			$pointer=1;
			print OUTPUT $_;
		}else{
			$pointer=0;
		}
	}else{
		if ($pointer==1) {
			print OUTPUT $_;
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
Usage:    search_fasta.pl -i inpute_file -d list -o output_file
Function: search a fasta file, and output all sequence ID in database
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -c column select, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-12-29
Notes:    Output input unique sequence
\n/
    )
}