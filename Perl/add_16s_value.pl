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
my %database; #database in hash
my %parent; #database in hash
while (<DATABASE>) {
	chomp;
#568	Root1480D1
#3	Root630_Root65
#1	Root136
	my @tmp=split/\t/;
	$database{$tmp[1]}=$tmp[0];
	$parent{$tmp[1]}=$tmp[1];
	print "$tmp[1]\t$tmp[0]\n";
	if ($tmp[1]=~/_/) {
		my @tmp1=split('_',$tmp[1]);
		foreach (@tmp1) {
			$database{$_}=$tmp[0];
			$parent{$_}=$tmp[1];
			print "$_\t$tmp[0]\n";
		}
	}
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[0]=' ' unless defined($tmp[0]);
	$database{$tmp[0]}=' ' unless defined($database{$tmp[0]});
	$parent{$tmp[0]}=' ' unless defined($parent{$tmp[0]});
	print OUTPUT "$tmp[0]\t$database{$tmp[0]}\t$parent{$tmp[0]}\n";
	
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
Usage:    add_16S_value.pl -i bac_OTU_list -o redup_bac_reads -d database	-h header num
Function: add vaule of each bacterial
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}