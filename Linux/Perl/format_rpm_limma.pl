#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my @list;
while (<DATABASE>) {
	chomp;
	next unless ($_=~/\w+/);
	my @tmp=split/\t/;
	push @list,$tmp[0];
}
close DATABASE;
my %gene;
foreach  $file(@list) {
	open DATABASE,"<$opts{i}$file";
	while (<DATABASE>) {
		chomp;
		my @tmp=split/\t/;
		$gene{$tmp[0]}{$file}=log($tmp[1]+1)/log(2);
	}
	close DATABASE;

}

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
foreach  (@list) {#output header
	print OUTPUT "\t$_";
}
print OUTPUT "\n";
#$line=0;
foreach $gi(sort keys %gene) {
#	$line++;
#	print "\n$line\n";
	print OUTPUT $gi;
	foreach $file(@list) {
		print OUTPUT "\t$gene{$gi}{$file}";
	}
	print OUTPUT "\n";
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
Usage:    format_rpm_limma.pl -i input_directory -d list -o output_file
Function: get each samples expression for limma analysis
Command:  -i inpute directory
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-04-02
Notes:    
\n/
    )
}