#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{o}="sRNA_basic.txt" unless defined($opts{o});


###############################################################################
#Main text.
###############################################################################
my %len;
my ($nr,$reads);
my %hash; #record whether miRNA exist
open INPUT,"<$opts{i}";
open OUTPUT,">>$opts{o}";
#if ($opts{t} eq "tab") {
	while (<INPUT>) {
		chomp;
		my @tmp=split/\t/;
		my $length=length($tmp[0]);
		$len{$length}{"nr"}++;
		$len{$length}{"reads"}+=$tmp[1];
		$nr++;
		$reads+=$tmp[1];
	}


print OUTPUT "length\tnonredundant\treads\n";
foreach (sort keys %len) {
	print OUTPUT "$_\t$len{$_}{nr}\t$len{$_}{reads}\n";
}
print OUTPUT "total\t$nr\t$reads\n";
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
Usage:    sRNA_basic.pl -i inpute_file -o sRNA_basic.txt
Function: Calculate sRNA length distribution, input tab or fasta
Command:  -i inpute file name (Must)
          -o output file name (Must), default sRNA_basic.txt
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-09-19
Notes:    
\n/
    )
}