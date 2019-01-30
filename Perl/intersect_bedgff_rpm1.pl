#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
my ($total,%reads,%rpm);
$total=$opts{t};
open DATABASE,"<$opts{d}";
#chrom0 Start1 End2			name3						score4 strand5	Chr6	source7				type8		start9	end10	score11	strand12 13		detail14(id,name,biotype)	overlap_length15
#Chr1    18090   18113   AGCAGATTACCATCACGATGGATG        2       +       Chr1    omicslab        SRC     17992   18167   5.52    +       .       ID=ath-SRC2;length_type=24-nt   23
#Chr1    14      37      ATTCAGAGGTTTAGGGTTTAGGGT        2       -       Chr1    omicslab        sRNA_cluster    8       1336    21.72   -       .       ID=SRC1;length_type=24-nt       23
#tur_scaffold192233      2575    2599    866456_1        255     +       tur_scaffold192233      omicslab        SRC     2599    2777    4.56    +       .       ID=AADD-SRC4603;length_type=24-nt       1
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[14]=~/ID=([^;]+)/;
	$src=$1;
	$tmp[3]=~/\w+\_(\w+)/;
	$reads{$src}+=$1;
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#1       omicslab        sRNA_cluster    3758    3861    6.57    -       .       ID=SRC1;length_type=24-nt
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+)/;
	if (defined($reads{$1})) {
		$rpm{$1}=$reads{$1}*1000000/$total;
		printf OUTPUT "$1\t%.2f\t$reads{$1}\n",$rpm{$1};
	}else{
		print OUTPUT "$1\t0\t0\n";
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
Usage:    intersect_bedgff_rpm1.pl -i cluster -d intersect_gff_result -o RPM -t total reads
Function: get the intersect result into RPM, reads, copies
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -t total reads
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
		  v2.0 add the reads count following the rpm
		  v3.0 add total reads, not need loci file, bed from bam file
Update:   2016-03-25
Notes:    
\n/
    )
}
