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
`mkdir $opts{o}`;
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
open OUTPUTL,">$opts{o}log\.txt";
while (<INPUT>) {
	if ($_=~/>/) {
		$_=~/>(\w+)/;
		open OUTPUT,">$opts{o}$1\.fa";
		print OUTPUTL "$opts{o}$1\.fa\n"
	}
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
Usage:    split_chr.pl -i fasta -o output_dir\/
Function: split fasta genome into each file, also create output_log file
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-08-22
Notes:    
\n/
    )
}