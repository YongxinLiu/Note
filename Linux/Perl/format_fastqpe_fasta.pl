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

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#@SRR488770.1 GSM918104_r1:5:1:564:495 length=36
#GATGACTGGGCACACATGCTGGCATCGTATGCCGTC
#+SRR488770.1 GSM918104_r1:5:1:564:495 length=36
#IIIIIIIIIIIIIIIIIIIIIIIIII?I3I(9I/=1
open OUTPUT,">$opts{o}";
open INPUT,"<$opts{i}";
<INPUT>;
while (<INPUT>) {
	push @tmp1,$_;
	<INPUT>;
	<INPUT>;
	<INPUT>;
}
close INPUT;
open INPUT,"<$opts{d}";
<INPUT>;
while (<INPUT>) {
	push @tmp2,$_;
	<INPUT>;
	<INPUT>;
	<INPUT>;
}
close INPUT;

for (0..$#tmp1) {
	print OUTPUT ">$_","f\n$tmp1[$_]>$_","r\n$tmp2[$_]";
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
Usage:    format_fastqpe_fasta.pl -i fastq -o fasta
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-05-04
Notes:    
\n/
    )
}
