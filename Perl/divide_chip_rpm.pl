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
$opts{n}=1 unless defined($opts{n});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
#Chr0	start1	end2	value3
#IGS-A	0	1	13.2266
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$name=$tmp[0].$tmp[1];
	$database{$name}=$tmp[3];
	#print $name,"\n";
}

close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#IGS-A	0	1	13.2266
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$name=$tmp[0].$tmp[1];
	#print $name,"\n";
	$database{$name}=$last if !defined($database{$name});
	$rpm=($tmp[3]+$opts{n})/($database{$name}+$opts{n});
	print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$rpm\n";
	$last=$database{$name};
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
Usage:    divide_chip_rpm.pl -i chip_file -o output_file -d input_file
Function: divide RPM
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
          -n norm value, default 1
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-08-03
Notes:    
\n/
    )
}