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
getopts( 'i:o:d:p:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{p}="bin_count_BS1.pl" unless defined($opts{p});


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr	start	end		ID	
#1	131665000 	135154000 	Cen1
#2	92910000 	95590000 	Cen2
#3	99794000 	100998000 	Cen3
#4	105356000 	106653000 	Cen4
#5	101978000 	106563000 	Cen5
open OUTPUT,">>$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT "$tmp[3]\t"; 
	print "$opts{p} -i $tmp[0]\_$opts{d} -o $opts{o} -c $tmp[0] -s $tmp[1] -e $tmp[2]\n";
	`$opts{p} -i $tmp[0]\_$opts{d} -o $opts{o} -c $tmp[0] -s $tmp[1] -e $tmp[2]`;
}


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
Usage:    batch_count_BS.pl -i position_list -d database_fixed -o mC_level -p bin_count_BS1.pl
Function: batch call bin_count_BS1.pl to caluclate some region mC level, input include chr, start and end, database each chr CGmap
Command:  -i inpute_file_list (Must)
          -o output_file
          -d CGmap of each chr
          -p call program
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-21
Notes:    
\n/
    )
}
