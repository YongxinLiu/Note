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
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#0chr	1start	2end		3				4			5			6	7end		8type	9mC	10strand	11
#5       2846    3246    Nucleosome:13   0.103415        5       3241    3242    CHH     0.1     G       1
#5       2846    3246    Nucleosome:13   0.103415        5       3243    3244    CHH     0.142857142857  C       1
#5       3041    3441    Nucleosome:14   0.103415        5       3042    3043    CHH     0.0     G       1
open OUTPUT,">$opts{o}";
#position	C	CG	CHG	CHH
#1	0.90	0.93	0.85	0.02
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my (%mC,%count); # bin total mC level
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$pos=$tmp[7]-$tmp[1];
	$mC{$pos}{C}+=$tmp[9];
	$count{$pos}{C}++;
	$mC{$pos}{$tmp[8]}+=$tmp[9];
	$count{$pos}{$tmp[8]}++;
}
foreach  $pos(sort {$a<=>$b;} keys %mC) {
	printf OUTPUT "$pos\t%.4f\t%.4f\t%.4f\t%.4f\n",$mC{$pos}{C}/$count{$pos}{C},$mC{$pos}{CG}/$count{$pos}{CG},$mC{$pos}{CHG}/$count{$pos}{CHG},$mC{$pos}{CHH}/$count{$pos}{CHH};
}


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
Usage:    nucleosome_point_BS.pl -i CGmap -o mC_average -c chromosome -s start -e end
Function: count each point 5mC level, position=5mC_end-nucl_start, average methyl level, include C, CG, CHG, CHH
Command:  -i inpute file name, *.CGmap (Must)
          -o output file name (Must)
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-06-04
Notes:    Output C, CG, CHG, CHH level
\n/
    )
}
