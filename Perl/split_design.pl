#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:c:p:', \%opts );
&usage unless ( exists $opts{i});
my $start_time=time;
$opts{o}=dirname($opts{i}) unless defined($opts{o});
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{c}=0 unless defined($opts{c});
$opts{p}="" unless defined($opts{p});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}

open OUTPUTL,">$opts{o}/list.txt";
while (<INPUT>) {
	chomp;
	if (/^([^\s]+)\s{0,}$/) {
		print $1,"\n";
		open OUTPUT,">$opts{o}/$opts{p}${1}.txt";
		print OUTPUTL $1,"\n";
	}elsif (/^\s{0,}$/) {
		next;
	}
	else{
		print OUTPUT $_,"\n";
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
        qq!
Usage:    split_chr.pl -i inpute_file -o output_file -c chrcolumn
Function: split genome into each chr one file in output directory
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -h header number, default 0
          -c column, set the chromsome site, default 0
          -p prefix
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.2
Update:   2017/12/25
Notes:    
\n!
    )
}