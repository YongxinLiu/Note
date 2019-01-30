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

###############################################################################
#Main text.
###############################################################################
my %len;
my ($nr,$reads);
open INPUT,"<$opts{i}";
open OUTPUT,">>$opts{o}";
while (<INPUT>) {
	chomp;
	next if />/;
	my @tmp=split/\t/;
	my $length=length($tmp[0]);
	$len{$length}{"nr"}++;
	$len{$length}{"reads"}+=$tmp[1];
	$nr++;
	$reads+=$tmp[1];
}
print OUTPUT "$opts{i}\nLength\tTotal\tUnique\tTotal%\tUnique%\n";
foreach (sort keys %len) {
	my $uniq_seq_per=$len{$_}{"nr"}/$nr*100;
	my $reads_per=$len{$_}{"reads"}/$reads*100;
	printf OUTPUT "$_\t$len{$_}{reads}\t$len{$_}{nr}\t%.2f\t%.2f\n",$reads_per,$uniq_seq_per;
}
print OUTPUT "Sum\t$reads\t$nr\n\n";
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
Usage:    sRNA_len.pl -i reads and count in tsf -o output_file
Function: Summary sRNA total reads and unique reads percentage in each length
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-12-05
Notes:    
\n/
    )
}
