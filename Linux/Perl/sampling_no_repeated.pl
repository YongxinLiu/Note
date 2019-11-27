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
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in array
my (@database);
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @database,$tmp[0];
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
my %sampling;
my $total=$#database;
while ((keys %sampling)<$opts{i}) { #random sampling, remove repeat by hash, until to set number
	$random=int(rand($total+1));
	$sampling{$database[$random]}=0;
}

foreach  (keys %sampling) {
	print OUTPUT "$_\n";
}

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
Usage:    sampling_no_repeated.pl -i sampling_number -o output_file -d database	-h header num
Function: Template for Perl
Command:  -i sampling_number (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-08-22
Notes:    
\n/
    )
}