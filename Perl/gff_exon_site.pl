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

open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";

my ($gene_id,$chromosome,$chain);
while (<INPUT>) {
	chomp;

###############################################################################
#Input exon and intron site
###############################################################################
	if (/gene/ && /Note/) {#get the gene chain and chromosome
		my @line=split('\t',$_);#0=chromosome,3=start_site,4=end_site,6=chain,8=gene_name;
		$chromosome=$line[0];
		$chain=$line[6];
	}
	if (/mRNA/) {#get one isoform start
		my @line=split('\t',$_);#0=chromosome,3=start_site,4=end_site,6=chain,8=gene_name;
		my @line8=split(';',$line[8]);
		$line8[0]=~s/ID=//;
		print OUTPUT "\n$chromosome\t$chain\t$line8[0]";
		$gene_id=$line8[0];
	}
	if (/exon/) {#"+"site from little to big, "-" site from big to little
		my @line=split('\t',$_);#0=chromosome,3=start_site,4=end_site,6=chain,8=gene_name;
		$line[8]=~s/Parent=//;
		if ($line[8] eq $gene_id){
			if ($chain eq "+"){
				print OUTPUT "\t$line[3]\t$line[4]";
			}else{
				print OUTPUT "\t$line[4]\t$line[3]";
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
