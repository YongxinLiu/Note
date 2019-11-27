#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:r:m:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /mnt/bai/yongxin/other/yangweicai/ath/temp/ColCvi24h1.bam
#1       4426    .       T       <*>     0       .       DP=1;I16=0,1,0,0,41,1681,0,0,60,3600,0,0,25,625,0,0;QS=1,0;MQ0F=0       PL:AD   0,3,41:1,0
#1       4648    .       C       <*>     0       .       DP=1;I16=1,0,0,0,32,1024,0,0,60,3600,0,0,25,625,0,0;QS=1,0;MQ0F=0       PL
open OUTPUT,">$opts{o}";
#chr0   start1  end2	3			4
#Chr1    8       31		Col_value	Cvi_value
while (<INPUT>) {
	chomp;
	if (/^#/) {
		next;
	}
	my @tmp=split/\t/;
	next unless defined($tmp[9]);
	$tmp[9]=~/:(\d+),(\d+)/;
#	print "$tmp[1]\t$1\t$2\n";
	next if $1+$2<1;
	print OUTPUT "$tmp[0]\t",$tmp[1]-1,"\t$tmp[1]\t$1\t$2\n";
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
        qq!
Usage:    format_vcf2bed.pl -i vcf_file -o bed_file
Function: format vcf file to bed file, for bedtools analysis
Command:  -i vcf_file
          -o bed file
          -d database file name (Opt)
#          -s scale, total mapped reads, normalize to RPM, default = 1 M
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/5/28
Notes:    
\n!
    )
}