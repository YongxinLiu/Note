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
open LOCI,"<$opts{d}";
#sequences0		reads1	loci2
#AAAAAAAAAAAAAACTCGGCC   1       1
my $total;
my %loci;
while (<LOCI>) {
	chomp;
	my @tmp=split/\t/;
	$total+=$tmp[1];
	$loci{$tmp[0]}=$tmp[2];
}
close LOCI;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
my (%reads,%reads_nrom);
while (<INPUT>) {
#chrom0 Start1 End2				name3			reads4 strand5 chrom6	Start7		End8 name9  score10 strand11 overlap12
#1       141012  141032  TATGTGTGTATCTTTCGGAAC   2       -       1       140907  141092  1       .       -       20
	chomp;
	my @tmp=split/\t/;
	$reads{$tmp[9]}+=$tmp[4];
	$loci{$tmp[3]}=1 unless defined($loci{$tmp[3]});
	$reads_norm{$tmp[9]}+=$tmp[4]/$loci{$tmp[3]};
}
close INPUT;
open OUTPUT,">$opts{o}";
#ID0		reads1	RPM2	reads_norm3	RPM_norm4		
#zma-miR1432-5p  65      7.49	65      7.49	
foreach  (sort keys %reads) {
		$rpm=$reads{$_}/$total*1000000;
		$rpm_norm=$reads_norm{$_}/$total*1000000;
		printf OUTPUT "$_\t$reads{$_}\t%.2f\t$reads_norm{$_}\t%.2f\n",$rpm,$rpm_norm;
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
Usage:    intersect_bed_rpm3.pl -d gff_cluster -i intersect_gff_result -t loci_database -o result
Function: use intersect result get assigned region RPM, not normalized and normalized, need genome total loci get total reads
Command:  -i inersect_bed result
          -o cluster or miRNA RPM
          -d  filtered_sRNA loci
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-08-30
Notes:    
\n/
    )
}