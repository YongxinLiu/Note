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
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;

#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#0		1		2		3start	4end	5	6 7 8ID
#Bd1	rep14	repeat	1	249	272	+	.	ID=bdi2-rep1;Note=Simple_repeat;Annotation=(CTAAACC)n
#Bd1	rep14	repeat	4487	4531	13	+	.	ID=bdi2-rep2;Note=Low_complexity;Annotation=A-rich
#Bd1	rep14	repeat	9736	9806	43	+	.	ID=bdi2-rep3;Note=Simple_repeat;Annotation=(AACCCTA)n
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my %repeat;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/Annotation=([^;]+)/;
	$repeat{$1}{num}++;
	$repeat{$1}{len}+=$tmp[4]-$tmp[3]+1;
}
close INPUT;

foreach  (keys %repeat) {
	print OUTPUT "$_\t$repeat{$_}{len}\t$repeat{$_}{num}\n";
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
Usage:    count_repeat_gff.pl -i repeat.gff -o repeat_length_number
Function: count repeatmasker gff result repeat type, length number, gff is 1 based coordinate
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-12-7
Notes:    
\n/
    )
}