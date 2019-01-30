#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:u:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{u}='F' unless defined($opts{u});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#10:25697367-25697548	1	89.67	184	16	1	1	181	37212238	37212055	5e-59	  231
#10:37558399-37558518	10	94.96	119	5	1	1	119	57519603	57519486	2e-45	  185
#10:37558004-37558097	10	94.68	94	4	1	1	93	57520687	57520594	3e-33	  145
open OUTPUT,">$opts{o}";
#chrom Start End name score strand; socre was replace by RPKM
#chr22 1000 5000 cloneA(length+identity) 960(score) + 
my $i=0;
my %seqid;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
#	if ($opts{u} eq 'T') {
#		next if defined($seqid{$tmp[0]});
#	}
	$i++;
	$tmp[0]=~/(\d+)\:(\d+)\-(\d+)/;
	$chr=$1;
	$start=$2;
	$end=$3;
	if ($tmp[9]>$tmp[8]) {# + strand
		print OUTPUT "Sbi$chr\t$start\t$end\t$tmp[1]\t$tmp[8]\t$tmp[9]\n";
	}else{
		print OUTPUT "Sbi$chr\t$start\t$end\t$tmp[1]\t$tmp[9]\t$tmp[8]\n";
	}
	$seqid{$tmp[0]}=0;
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
Usage:    format_blastn6_bed_pair.pl -i blastn table format -o bed
Function: format blastn table format to bed, ID is the query name ID, only keep chr ID
Command:  -i inpute blastn file name (Must)
          -o output bed file name (Must)
          -d database file name
		  -u unique, only report best 1 blast site, default=F, can set T
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-07-19
Notes:    2014-08-28	Add -u unique result 
\n/
    )
}