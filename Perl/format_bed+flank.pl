#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:f:r:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my %database;
#ID		length
#2000.1  243767
#2000.10 16101
#2000.11 6505
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[1];
}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr-	start1	end2	ID3	socre4	strand5
#tur_scaffold93179	11154	11254	1_9033	  171	+
#ata_Scaffold8258	1808	1895	2_8779	  158	+
#ata_Scaffold20541	2626	2725	3_7624	  185	+
open OUTPUT,">$opts{o}";
#chr-	start1	end2	ID3	socre4	strand5
#tur_scaffold93179	11154	11254	1_9033	  171	+
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$start=$tmp[1]-$opts{f};
	$start=0 if $start<0;
	$end=$tmp[2]+$opts{f};
	$end=$database{$tmp[0]} if $end>$database{$tmp[0]};
	print OUTPUT "$tmp[0]\t$start\t$end\t$tmp[3]\t$tmp[4]\t$tmp[5]\n";

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
Usage:    format_bed+flank.pl -i bed(line6) -o bed -d database_length
Function: add flank to bed file
Command:  -i input
          -o output file name (Must)
          -d database file name (Must)
          -f flank
          -r remove duplication
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-11-02
Notes:    
\n/
    )
}