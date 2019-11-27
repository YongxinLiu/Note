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
#database in hash
my %database; 
while (<DATABASE>) {#read each database total count
#ID	total	filter	percentage
#GSM424847       1091917 304739  27.91%
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[2];
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
my @filelist=glob "$opts{i}";
foreach $file(@filelist){
	open INPUT,"<$file";
	$file=basename($file);
	open OUTPUT,">$opts{o}$file";
	$file=~/(^[\d\w\_]+)/;
	print "$1\n";
	$total=$database{$1};
	while (<INPUT>) {
#chr	start	end	value
#Chr1    10      10      -9
		my @tmp=split/\t/;
		$value=$tmp[3]*1000000/$total;
		printf OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t%.2f\n",$value;
	}
	close DATABASE;
}
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	
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
Usage:    normalize.pl -i inpute_file -o output_file -d database	-h header num
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-08-22
Notes:    
\n/
    )
}