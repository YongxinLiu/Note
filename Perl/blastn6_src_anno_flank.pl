#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:q:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";#reads fasta in hash
#Chr1    omicslab        sRNA_cluster    17992   18166   36.55   +       .       ID=SRC2;length_type=24-nt       AT1G01030|Chr1|11649|13714|-|4277|TAIR10|protein_coding_gene|AP2/B3-like transcriptional factor family protein;ath-rep20|Chr1|17010|17256|+|735|Repbase2012|RC/Helitron|ATREP3  AT1G01040|Chr1|23146|31227|+|4979|TAIR10|protein_coding_gene|dicer-like 1;ath-rep22|Chr1|18661|18728|+|494|Repbase2012|RC/Helitron|ATREP20
my %database; #open ath database fasta
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
#	print "$tmp[8]\t$tmp[9]\t$tmp[10]\n";
	$database{$tmp[8]}{left}=$tmp[9];
	$database{$tmp[8]}{right}=$tmp[10];
}
close DATABASE;
my %query; #open osa query fasta
open DATABASE,"<$opts{q}";#reads fasta in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$query{$tmp[8]}{left}=$tmp[9];
	$query{$tmp[8]}{right}=$tmp[10];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#IDq0	IDd1	similiraty2	length3	mismatch4	gap5	qstart6	qend7	dstart8	dend9	E-value10	score11
#IDq	IDd			92.59   27			 2       0       296     322     885     911     0.004   39.9
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT "$_\t$query{$tmp[0]}{left}\t$query{$tmp[0]}{right}\t$database{$tmp[1]}{left}\t$database{$tmp[1]}{right}\n";
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
Usage:    blastn6_src_anno_flank.pl -i blastn6format -o output.tsv -d database.flank -q query.flank
Function: add blast result query and target flank anno (query left, query right, target left, target right)
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file fasta
          -q query file fasta
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-04-15
Notes:    
\n/
    )
}