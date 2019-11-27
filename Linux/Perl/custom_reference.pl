#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:f:g:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{f}=100000 unless defined($opts{f});


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}\.bed";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
open TMP,">tmp.bed";
my $start;
my $end=(-1)*$opts{f};
while (<INPUT>) {
#chr start	end id
#3	218920000	219270000	Dp3a	.	+
#9	55457000	56180000	sDic15	.	+
	chomp;
	my @tmp=split/\t/;
	print TMP "$tmp[0]\t",$tmp[1]-$opts{f},"\t",$tmp[2]+$opts{f},"\t$tmp[3]\n";
	$len=$tmp[2]-$tmp[1];
	$start=$end+2*$opts{f};
	$end=$start+$len;
	print OUTPUT "$opts{o}\t$start\t$end\t$tmp[3]\n";
	}
close INPUT;
close OUTPUT;
`fastaFromBed -fi $opts{g} -bed tmp.bed -fo tmp.txt -name -tab`;

open INPUT,"<tmp.txt";
my $seq;
open OUTPUT,">$opts{o}\.fa";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$seq.=$tmp[1];
}
print OUTPUT ">$opts{o}\n$seq\n";
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
Usage:    custom_reference.pl -i inpute_file.bed -o output.fa&bed -f flank -g genome file
Function: based on a bed file create reference fasta and bed. first add flank on bed, then get fasta from genome, merge fasta as new chromosome, create bed file for new chromosome
Command:  -i bed (Must)
          -o fa and bed (Must)
          -f flank region default 100kb
          -g genome file
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-09-27
Notes:    
\n/
    )
}