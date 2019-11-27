#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:l:r:g:q:t:', \%opts );
&usage unless ( exists $opts{l} && exists $opts{r} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{l} and $opts{r}\nOutput file is $opts{o}\.bam\n";
print "Genome file is $opts{g}\n\n" if defined($opts{g});
$opts{t}=10 unless defined($opts{t});
$opts{q}=20 unless defined($opts{q});

print "Mapping left read to reference\n";
print "bwa aln -t $opts{t} $opts{g} $opts{l} -f $opts{o}\_left.sai\n\n";
`bwa aln -t $opts{t} $opts{g} $opts{l} -f $opts{o}\_left.sai`;

print "Mapping right read to reference\n";
print "bwa aln -t $opts{t} $opts{g} $opts{r} -f $opts{o}\_right.sai\n\n";
`bwa aln -t $opts{t} $opts{g} $opts{r} -f $opts{o}\_right.sai`;

print "Merge pair-end mapping and output sam\n";
print "bwa sampe $opts{g} $opts{o}\_left.sai $opts{o}\_right.sai $opts{l} $opts{r} -f $opts{o}\_pe.sam -n 1 -N 1\n\n";
`bwa sampe $opts{g} $opts{o}\_left.sai $opts{o}\_right.sai $opts{l} $opts{r} -f $opts{o}\_pe.sam -n 1 -N 1`; #-n 1 -N 1 for exclude repeat, only for chip-seq

#print "Delete temp sai file\n\n";
#`rm $opts{o}\_left.sai`; 
#`rm $opts{o}\_right.sai`; 

print "Statistic mapping result\n";
print "samtools view -bS $opts{o}\_pe.sam|samtools flagstat - \n\n";
`samtools view -bS $opts{o}\_pe.sam|samtools flagstat - >$opts{o}\.log 2>&1`; #no statistic on screen

print "Sam to sorted bam, filter unmapped, quality <20 result\n";
print "samtools view $opts{o}\_pe.sam -F4 -bS -q 20|samtools sort -  $opts{o}\n\n";
`samtools view $opts{o}\_pe.sam -F4 -bS -q 20|samtools sort -  $opts{o}`;

print "Index bam\n";
print "samtools index $opts{o}\.bam\n\n";
`samtools index $opts{o}\.bam`;


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
Usage:    mapping_bwa_pe.pl -l left.fq -r right.fq -g genome -t thread -q mapping_quality -o output_file_prefix
Function: mapping clean pair-end fastq to genome by bwa
Command:  -l leaf.fa (Must)
          -r right.fa (Must)
          -g reference genome
          -o output filtered unmapped, sorted bam and index
          -t threads, default 10
          -q mapping quality, default 20
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-06-06
Notes:    
\n/
    )
}