#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{t}=1.3 unless defined($opts{h});

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
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#0SRCs_ID	1length_type	2distrance_with_gene,less<1mean_in_gene 3gene_ID	4salt_leaf_rpm	5normal_leaf_rpm	6fold	7normal_leaf_FPKM	8salt_leaf_FPKM
#ath-SRC16	24-nt	-1535	AT1G01070	1.82	3.81	0.59	7.78383	10.006
	chomp;
	my @tmp=split/\t/;
	next if ($tmp[6]>1/$opts{t} && $tmp[6]<$opts{t}); #filter SRCs not changed
	next if !defined($tmp[8]); #filter related gene not expreseed
	$fold=($tmp[8]+1)/($tmp[7]+1);
	next if ($fold>1/$opts{t} && $fold<$opts{t});#filter gene not changed
	if ($tmp[6]>$opts{t} && $fold>$opts{t}) {
		printf OUTPUT "$_\t%.2f\tpositive\n",$fold;
	}elsif ($tmp[6]<1/$opts{t} && $fold<1/$opts{t}) {
		printf OUTPUT "$_\t%.2f\tpositive\n",$fold;
	}elsif ($tmp[6]>$opts{t} && $fold<1/$opts{t}) {
		printf OUTPUT "$_\t%.2f\tnegative\n",$fold;
	}elsif ($tmp[6]<1/$opts{t} && $fold>$opts{t}) {
		printf OUTPUT "$_\t%.2f\tnegative\n",$fold;
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
Usage:    SRCs_mRNA_relation.pl -i inpute_file -o output_file -d database	-h header num
Function: Analysis SRCs and related gene expression, and find positive or negtive regulate SRCs
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
          -t threshold of fold change, default 1.3
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-08-22
Notes:    
\n/
    )
}