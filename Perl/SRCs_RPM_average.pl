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
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

##database in array
#my (@tmp1,@tmp2);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;
open INPUT,"<$opts{d}";
my @filelist;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print $tmp[0],"\n";
	push @filelist,$tmp[0];
	
}
close INPUT;

my %list;
foreach $file(@filelist){
	open DATABASE,"<$opts{i}$file";
	$file=basename($file);
	while (<DATABASE>) {
		my @tmp=split/\t/;
		$list{$tmp[0]}{$file}=$tmp[1];
	}
	close DATABASE;
}

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
foreach  (sort keys %list) {
	$average=&average($_);
	printf OUTPUT "$_\t%.2f\n",$average;
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
Usage:    SRCs_RPM_average.pl -i input_dir	-d list -o output_file
Function: get a list samples SRCs RPM average 
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-22
Notes:    
\n/
    )
}
sub average{
	my $tmp=$_[0];
	my $sum;
	@key2=sort keys %{$list{$tmp}};
	foreach  (@key2) {
		$sum+=$list{$tmp}{$_};
	}
	$average=$sum/($#key2+1);
}