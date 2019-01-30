#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Bio::SeqIO;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################

###############################################################################
#Read the genome in hash
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database;
#my $chr;
#while (<DATABASE>) {
#	chomp;
#	if (/>/) {
#		$chr=(split('\s',$_))[0];
#		$chr=~s/>//;
##		print $chr,"\n";
#	}else{
#		$database{$chr}.=$_;
#	}
#}
#close DATABASE;
##print substr($database{$chr},0,100),"\n",substr($database{$chr},-100,100);

open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
#open OUTPUTS,">>$opts{s}";
my ($gene_id,$start_site,$end_site,$chain);
while (<INPUT>) {
	chomp;


###############################################################################
#Distill transcript from genome
###############################################################################
#	if (/gene/ && /Note/) {
#		#Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
#		my @line=split('\t',$_);#0=chromosome,3=start_site,4=end_site,6=chain,7=gene_name;
##		print $line[0],"\n",$line[3],"\n",$line[4],"\n",$line[6],"\n",$line[8],"\n";
#		#ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
#		my @line8=split(';',$line[8]);
#		$line8[0]=~s/ID=//;
#		print OUTPUT ">$line8[0]\t$line8[1]\n";
##		print $line8[0],"\n";
#		$line[0]=~s/Chr//;
#		foreach  (keys %database) {
#			if ($_=~m/^$line[0]/i) {
#				my $sequence=substr($database{$_},$line[3]-1,$line[4]-$line[3]+1);
#				if ($line[6] eq "+") {
#					print OUTPUT $sequence,"\n";
#				}elsif ($line[6] eq "-"){
#					$sequence=reverse $sequence;
#					$sequence=~tr/ACGTacgt/TGCAtgca/;
#					print OUTPUT $sequence,"\n";
#				}
#			}
#		}
#	}

	
###############################################################################
#Input exon and intron site
###############################################################################
	if (/gene/ && /Note/) {
		my @line=split('\t',$_);#0=chromosome,3=start_site,4=end_site,6=chain,8=gene_name;
		$start_site=$line[3];
		$end_site=$line[4];
		$chain=$line[6];
	}
	if (/mRNA/) {
		my @line=split('\t',$_);#0=chromosome,3=start_site,4=end_site,6=chain,8=gene_name;
		my @line8=split(';',$line[8]);
		$line8[0]=~s/ID=//;
		print OUTPUT "\n$line8[0]";
		$gene_id=$line8[0];
	}
	if (/exon/) {
		my @line=split('\t',$_);#0=chromosome,3=start_site,4=end_site,6=chain,8=gene_name;
		$line[8]=~s/Parent=//;
		if ($line[8] eq $gene_id) {
			if ($chain eq "+") {
				print OUTPUT "\t",$line[3]-$start_site+1,"\t",$line[4]-$start_site+1;
			}else{
				print OUTPUT "\t",$end_site-$line[4]+1,"\t",$end_site-$line[3]+1;
			}
		}
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
Usage:    template.pl -i inpute_file -o output_file
Function:
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-09-20
Notes:    
\n/
    )
}