#!/usr/bin/perl
use warnings;
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{t}="int" unless defined($opts{t});
###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split('\t',$_);
#}
#close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>){
	chomp;
	my @tmp=split('\t',$_);
	print OUTPUT "$tmp[0] char(30) not null primary key,";
	for (1..$#tmp){
		print $_;
		print OUTPUT "$tmp[$_] $opts{t},";
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
Usage:    mysql_title.pl -i inpute_file -o output_file -t database type
Function: format title file to mysql used title format
Command:  -i inpute file name (Must)
          -o output file name
          -t database type, default int
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-16
Notes:    
\n/
    )
}
