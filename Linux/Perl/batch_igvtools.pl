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
getopts( 'i:o:t:g:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{t}="toTDF" unless defined($opts{t});
$opts{g}="Zma3" unless defined($opts{g});


###############################################################################
#Main text.
###############################################################################
my @filelist=glob "$opts{i}";
for (0..$#filelist){
	my $basename=basename($filelist[$_]);
	my $output=$opts{o}.$basename."\.tdf";
	print "~/software/IGVTools/igvtools $opts{t} $filelist[$_] $output $opts{g}\n";
	`~/software/IGVTools/igvtools $opts{t} $filelist[$_] $output $opts{g}`;
}
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
Usage:    batch_igvtools.pl -i inpute_file_list -o output_file_directory -t toTDF -g tair10
Function: batch call igvtools toTDF
Command:  -i inpute_file_list (Must)
          -o output_file_directory
          -t igvtools function type, such as toTDF
          -g genomeID, such as tair10, hg18, Zma3
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-04
Notes:    IGVtools from PCCE-hdsu
\n/
    )
}
