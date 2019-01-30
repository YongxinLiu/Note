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

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#>618 TCONS_00000619 C166445605. 15-681
#CGTGCCCATGCCGCAGCTCCTCTAGGACCAGCCTCCATTGACAATAGATAGAGTCTAGGC
#TTGAGCACGGATGCAAACAGCAGCATGCCGTCCTCCACACCGATGAGGTAAAGGAATCCC
#TGGCCCAGGTTCTCTGGCACATTGATCACAGATAGTTGTTGCTCACGAATGTTGAACTCC
#ACGAGACTGTCACACTCATAACTGGGGACATAGACCTTATCTCCCACCACAGCAGTGTGC
#CCTCCGTCGTCAGTGTGCAAACCCGGGCTCCGAACAGAGATCTTGTCGCTCCACTGACAA
#GTCTCCGACGAGTAGACAGTGGCGAGGGCAATGCCTCCTTCCTCCACAGAGCCCACCAAG
#GCCACGAGGAAAGGACCGCCGCCATGGCAGTCGAGATGGCCACATCGATCTTTGGCGCAG
#AGCGCCGTGGCATTACAGCGTATCCCAGTCCCGCTGAAGAAGTCGTCGGGGCCGTCGCCC
#GGCCAGTGCATGATGTCGCGGCACTTGGGGTCGTCGTCGATGTCCCACCAGTTGCCGGTG
#ACGAGGTCGCGGGCCTCAAAAGGGGGGTCGTACGCAGGGGCGTAGAGGAGGACGAGGCCG
#TGCCGGGAGTCGAGGACCTGGCAGTCTGGGCGGTCCTGGGACGGCCGGCGGCGGAAGGTC
#GAAGTGG
#>619 TCONS_00000620 C166447525- 6-203,288-476
#GATGTCGCCCTTGCCCTTGCCGGCCTTCTTGGGGAGCAGCGTGGTGTGGATGCTGGGCAG
#CACCCCGCCGGCGGCGATGGTGACCGAGCCCAGCAGCCGGCTCAGCTCCTCGTCGTTGCG
#CACCGCCAGCTGGATGTGGCGCGGCACGATCCGGTTCTTCTTGTTGTCCCGCGCCGCGTT
#CCCGGCGAGCTCCAGGGTCTCTGCGGCGAGGTACTCGAGGACGGCGGCGAGGTAGACGGG
#GGCGCCGGCGCCGACGCGCTGCGCGTACTTGCCGACCTTGAGGTACCGGGCGACGCGGCC
#GACGGGGAACTGCAGGCCGGCCTTGGACGACCGGGACACGGACTTGGCGCTGGCCGCGGG
#CTTCGCCTTGCCCCTCCCACTGCCGGC
#>620 TCONS_00000621 C166447639+ 1-398,560-626
#GAGTGCCGCCGTGCCCTGCTCCCTCCTCTCTTCCTCTGTCTCCTCTCACCTCGCTCATCC
#CCCAATCTCTCTGTCCCCAGGACGAACAGGAGCATGCCATGGCCTCGGATCCGCGCCGCC
#GTCGCCGGAGAACGTCGCCCGCCACCGGTTCTTCGACCACAACTACCGGATCAAGTTGCC
#CCGCCTCCGATCTACCCGTTCCCGGCCACCCCAACCTCGCCTCGTCATTGGGCGCGCAGA
#GAGACGGCCATGGCCGCCCGATTACCGAGCGCCGCCGGCCAAGGTGGCTGCCCGATGCAC
#GAGCGCCGCCGGCCAAGGACCCGCCGGCGTCCAGATCCCCTCCACGCCATCCCCATCTCG
#CCGGAGACCTGCTCCAAGTCTGGCCGCCCTCCAAACAGTACGCGGAGACAGAGGCCAAGG
#ATGACCGGCAACTCATGGAGGAGGACACCTCAAATCAGTCTGGGA
open OUTPUT,">$opts{o}";
my %seq;
while (<INPUT>) {
	chomp;
	if (/>/) {
		$name=$_;
	}
	if (/(^[AGCTN]+$)/) {
		$seq{$name}.=$_;
	}
}
foreach  (keys %seq) {
	next if /\.\s/;
	if (/\+\s/) {
		print OUTPUT "$_\n$seq{$_}\n";
	}
	if (/\-\s/) {
		$seq=$seq{$_};
		$seq=~tr/ACGT/TGCA/;
		$seq=reverse($seq);
		print OUTPUT "$_\n$seq\n";
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
Usage:    format_fasta_revcom.pl -i inpute_file -o output_file
Function: gtf_to_fasta, reads on minus strand need to be revcom, this script for this function
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-06-06
Notes:    
\n/
    )
}