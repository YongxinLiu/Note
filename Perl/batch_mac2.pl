#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:s:q:l:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{h}=0 unless defined($opts{h});
$opts{i}="doc/IPIn.txt" unless defined($opts{i});
$opts{o}="result/" unless defined($opts{o});
$opts{s}=1000000 unless defined($opts{s});
$opts{q}=0.05 unless defined($opts{q});
$opts{l}="result/log.txt" unless defined($opts{l});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print "macs2 callpeak -t result/$tmp[0].bam -c result/$tmp[1].bam -f BAM -g $opts{s} --outdir $opts{o} --name $tmp[0] --bdg --qvalue $opts{q} 2>$opts{l}\n";
	`macs2 callpeak -t result/$tmp[0].bam -c result/$tmp[1].bam -f BAM -g $opts{s} --outdir $opts{o} --name $tmp[0] --bdg --qvalue $opts{q} 2>$opts{l}`;
}
close INPUT;

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
Usage:    macs2 callpeak -t ChIP.bam -c Input.bam -f BAM -g genome_size --outdir output --name ChIP --bdg --qvalue 0.95 2>log.txt
Function: batch call program have one input, and one output, and an database according with input
		  default include double strand, at least 1bp overlap
Command:  -i inpute_file, include samples list
          -o output_file_directory
          -s genome size
          -q qvalue, default 0.05
          -l log file, default result/log.txt
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/6/11
Notes:    
\n!
    )
}
