#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:p:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################
my @filelist=glob "$opts{i}";
#my @filelist=glob "~/data/Data/sRNA/GSM/*";
foreach (@filelist){
	my $basename=basename($_);
	my $output=$opts{o};
	print "$opts{p} -i $_ -o $output\n";
	`$opts{p} -f $_ -o $output`;
#	print $opts{o},"\n";
#	print $_,"\n";
}

#print $#filelist,"\n"
print scalar @filelist," files have been treated.\n";


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
Usage:    batch1.pl -i inpute_file_list -o output_file_directory -p call program
Function: batch call program have one input, and one output
Command:  -i inpute_file_list (Must)
          -o output_file_directory
          -d call program
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-06-03
Notes:    
\n/
    )
}
