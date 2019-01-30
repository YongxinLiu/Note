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
#qseqid0								sseqid1	pident2 length3 mismatch4 gapopen5 qstart6 qend7    sstart8	send9		evalue9 bitscore10
#RF00001;5S_rRNA;L42764.1/1575-1690      9       100.00  116     0	    0       1       116     131258636	131258751	3e-54    215
#0		1		2		3		chr4	start5	end6		7			strand8	description9	biotype10		11		12		13		14
#701     32.8    6.3     4.7     1       2121    2893    (301351242)     +       LINE1-2_ZM      LINE/L1 3034    4148    (317)   1
open OUTPUT,">$opts{o}";
#id0	chr1 strand2	start3	end4 	length5
#1       1       -       3758    3861    104   
my $i=0;
while (<INPUT>) {
	chomp;
	$i++;
	my @tmp=split/\t/;
	if ($tmp[9]>$tmp[8]) {# + strand
		print OUTPUT "$i\t$tmp[1]\t+\t$tmp[8]\t$tmp[9]\t",$tmp[9]-$tmp[8]+1,"\n";
	}else{
		print OUTPUT "$i\t$tmp[1]\t-\t$tmp[9]\t$tmp[8]\t",$tmp[8]-$tmp[9]+1,"\n";
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
Usage:    template.pl -i inpute_file -o output_file
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-05-04
Notes:    
\n/
    )
}