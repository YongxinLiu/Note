#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open LOCI,"<$opts{t}";
##sequences0		reads1	loci2
##AAAAAAAAAAAAAACTCGGCC   1       1
#my (%loci,$total);
#while (<LOCI>) {
#	chomp;
#	my @tmp=split/\t/;
#	$loci{$tmp[0]}=$tmp[2];#preserve all the loci number in hash
#	$total+=$tmp[1];
#}
#close LOCI;
#open DATABASE,"<$opts{d}";
##chrom0	source1	type2 start3  end4		5 strand6		7		ID8															
##9       ensembl gene    66347   68582   .       -       .       ID=GRMZM2G354611;Name=GRMZM2G354611;biotype=protein_coding
#my (%id,%len);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$tmp[8]=~/Name=([^;]+)/;
##	print $1,"\n";
#	$id{$1}=$1;
#	$len{$1}=$tmp[4]-$tmp[3]+1;
#}
#close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
my %repeat;
my %family;
while (<INPUT>) {
#chrom0 start1			end2		chrom3		source4	type5	start6			end7			8	 strand9	10		ID11												overlap12
#1      36628326        36628411        1       rep14   repeat  36627855        36630708        19840   +       .       ID=bdi-rep14284;Note=LTR/Gypsy;Annotation=Gypsy-18_BD-I 85
#1       36628412        36628413        1       rep14   repeat  36627855        36630708        19840   +       .       ID=bdi-rep14284;Note=LTR/Gypsy;Annotation=Gypsy-18_BD-I 1
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[11]=~/Note=([^;]+);Annotation=([^;]+)/;
	$repeat{$2}+=$tmp[12];
	$family{$1}+=$tmp[12];
}
close INPUT;
open OUTPUT,">$opts{o}";
foreach  (sort {$family{$b}<=>$family{$a};} keys %family) {
	print OUTPUT "$_\t$family{$_}\n";
}
print OUTPUT "\n\n\n";
foreach  (sort {$repeat{$b}<=>$repeat{$a};} keys %repeat) {
	print OUTPUT "$_\t$repeat{$_}\n";
}

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
Usage:    intersect_bedgff_repeat_stat.pl -i intersect_bedgff -o repeat_length
Function: Calculate each repeat family coverage
Command:  -i inersect_bed result
          -o each repeat family length
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-11-19
Notes:    
\n/
    )
}