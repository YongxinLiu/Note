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
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[0];
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
#0	1			2	3			4			5		6	7	8									9
#3	omicslab	SRC	218848248	218848418	1139.32	-	.	ID=zma-SRC2214482;length_type=other	LTR/opie|RLC_opie_AC197691_5727|3|218848417|218849555|+|2
	my @tmp=split/\t/;
	print OUTPUT "$_\n" if defined($database{$tmp[0]});
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
Usage:    filter_gff_chr.pl -i gff -o filtered_gff -d database	-h header num
Function: filter gff by chromosome ID
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