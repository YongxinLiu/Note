#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
#$opts{m}=-1 unless defined($opts{m});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr0	start1	end2	sequences3			reads4		strand5
#IGS-A   10      34      5998225_2       255     -
open OUTPUT,">$opts{o}"; 
my %count;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[3]=~/([\w\.]+)_/;
	$count{$tmp[0]}{$1}++;
}
# print title
@tmp=keys %count;
print OUTPUT "ID";
foreach my $key2 (sort keys %{$count{$tmp[0]}}) {
		printf OUTPUT "\t$key2";
}
print OUTPUT "\n";
@tmp=sort keys %{$count{$tmp[0]}};

foreach  my $key1 (sort keys %count) {
	print OUTPUT $key1;
	foreach my $key2 (@tmp) {
		$count{$key1}{$key2}=0 unless defined($count{$key1}{$key2});
#		$count{$key1}{$key2}=0 unless $count{$key1}{$key2}>0;
		print OUTPUT "\t$count{$key1}{$key2}";
	}
	print OUTPUT "\n";
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
Usage:    format_bed2otutable.pl -i bedfrombam -o otu_table
Function: result in one file, not normalize by loci, normalized to RPM
Command:  -i use filter loci bed
          -o output file name (Must)
          -d database file name (Opt)
          -a min length, 20 (Opt)
          -b min length, 25 (Opt)
          -s scale, total mapped reads, normalize to RPM, default = 1 M
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.2
Update:   2016-03-24
Notes:    Output chr by input
\n/
    )
}