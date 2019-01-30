#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:a:d:r:h:c:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{r}=1 unless defined($opts{r});
$opts{c}=9 unless defined($opts{r});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#tracking_id0	class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	AA1_FPKM9	AA1_conf_lo	AA1_conf_hi	AA1_status	AA2_FPKM13    AA2_conf_lo	AA2_conf_hi	AA2_status
#EPlTURT00000000001	-	-	EPlTURG00000000001	MIR1122	-	scaffold9027:9977-10070	93	-	0	0	0	OK	0	0	0	OK
open DATABASE,"<$opts{d}";
my %database; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[$opts{c}];
}
close DATABASE;
open DATABASE,"<$opts{a}";
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[$opts{c}];
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#IDA	IDD	score
#EPlTURT00000001365	EPlATAT00000000383	  114
#EPlTURT00000002503	EPlATAT00000001682	  157
#TRIUR3_00024-T1	EMT27083	  161
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	printf OUTPUT "$tmp[0]\t$tmp[1]\t$database{$tmp[0]}\t$database{$tmp[1]}\t%.4f\n",($database{$tmp[1]}+1)/($database{$tmp[0]}+1);
	
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
Usage:    cuffdiff_homolog1.pl -i homolog_list -o homolog_gene_exp -a AA_iso_exp -d DD_iso_exp -h header num
Function: Calculate homolog different, can select samples line 9, 13, 17
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -a databaseA, such as AA isoform expression from cuffdiff
          -d database file name, such as DD isoform expression from cuffdiff
          -r replication, default 1 means no rep, can set 2,3,4....
          -c column line, 9, 13, 17 represent sample1, 2, 3
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2014-09-16
Notes:    
\n/
    )
}