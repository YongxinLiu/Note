#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:p:d:e:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}"; # read homolog gene
#ID0	ID1	score2
#EPlTURT00000001365	EPlATAT00000000383	  114
#EPlTURT00000002503	EPlATAT00000001682	  157
my %homolog; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$homolog{$tmp[0]}=0;
	$homolog{$tmp[1]}=0;
}
close DATABASE;

open DATABASE,"<$opts{e}"; # read half homolog gene
my %homolog1; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$homolog1{$tmp[0]}=0 unless defined($homolog{$tmp[0]});
	$homolog1{$tmp[1]}=0 unless defined($homolog{$tmp[1]});
}
close DATABASE;

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
open OUTPUTP,">$opts{p}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) { # read transcripts, output half homolog and specific gene list
	chomp;
#chr0	start1	end2	ID3
#Scaffold20028	1816	5218	EMT23773	0	+	1816	5216	0   108,117,1079,219,37,	0,216,1531,2808,3366,
#Scaffold20028	19133	19258	EPlATAT00000000349	0	+	19133	1925126,	0,
	my @tmp=split/\t/;
	next if defined($homolog{$tmp[3]});
	if (defined($homolog1{$tmp[3]})) {
		print OUTPUT "$tmp[3]\n";
		next;
	}
	print OUTPUTP "$tmp[3]\n";
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
Usage:    AD_gene_half_unique.pl -i inpute_file -o output_file -d homolog gene -e half homolog gene	-h header num
Function: base on best_hit_mutual.pl result, find half homolog and specific genes of AADD
Command:  -i bed12 from gtf
          -o output multiple gene members ID
		  -p output unique gene ID
          -d homolog gene pair
          -e half homolog gene pair
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-09-21
Notes:    
\n/
    )
}