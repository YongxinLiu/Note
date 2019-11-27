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
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#qseqid0				sseqid1	pident2 length3 mismatch4 gapopen5 qstart6 qend7    sstart8	send9		evalue10 bitscore11
#IGS-A-long      scaffold51788   92.67   2700    173     23      122     2810    7463    10148   0.0     3866
#IGS-A-long      scaffold51788   91.41   1641    123     18      122     1756    17186   18814   0.0     2233
open OUTPUT,">$opts{o}";
#chrom Start End name score strand; socre was replace by RPKM
#chr22 1000 5000 cloneA(length+identity) 960(score) + 
my $i=0;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next if defined($hash{$tmp[0]});
	$hash{$tmp[0]}=$tmp[1];
	print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[10]\t$tmp[11]\n";
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
Usage:    format_blastn6_genepair.pl -i blastn table format -o bed
Function: base on blast_m6 result, and select each gene best target. output include: query seq, target seq, Evalue and bitscores
Command:  -i inpute blastn file name (Must)
          -o output bed file name (Must), include ID, target ID, Evalues, and scores
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-01
Notes:    
\n/
    )
}