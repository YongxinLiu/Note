#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#get the gff format annotation
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+);/;
	$database{$1}=$_;
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my %count;
my @list;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT "$database{$tmp[0]}\n";
	$count{miRNA}++ if $database{$tmp[0]}=~/miRNA/;
	$count{Repbase}++ if $database{$tmp[0]}=~/Repbase/;
	$count{protein_coding_gene}++ if $database{$tmp[0]}=~/protein_coding_gene/;
	$count{intergenic}++ if $database{$tmp[0]}=~/intergenic/;
	@tmp1=split('\t',$database{$tmp[0]});
	@tmp2=split('\;',$tmp1[9]);
	foreach  (@tmp2) {
		if (/protein_coding_gene/) {
			@tmp3=split/\|/;
			push @list,$tmp3[0];
		}
	}
}
foreach  (@list) {
	print OUTPUT "$_\n";
}

print OUTPUT "miRNA\t$count{miRNA}\nRepeat\t$count{Repbase}\nProtein\t$count{protein_coding_gene}\nintergenic\t$count{intergenic}\n";


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
Usage:    search_SRC_list.pl -i SRCs_list -d annotation -o SRCs_anno
Function: search a SRCs list for annotation
Command:  -i SRCs list in the file first line
          -o output SRCs annotation
          -d SRCs annotation file, usually gff format, contained SRCID with semicolon close
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-10
Notes:    
\n/
    )
}