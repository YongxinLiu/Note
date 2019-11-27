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
open LOCI,"<$opts{t}";
#sequences0		reads1	loci2
#AAAAAAAAAAAAAACTCGGCC   1       1
my (%loci,$total);
while (<LOCI>) {
	chomp;
	my @tmp=split/\t/;
	$loci{$tmp[0]}=$tmp[2];#preserve all the loci number in hash
	$total+=$tmp[1];
}
close LOCI;
open DATABASE,"<$opts{d}";
#id0	chr1 strand2	start3	end4 	length5	loci6	reads7	reads_nrom8(normalized,reads/loci)	rpm9(reads_norm/KM)
#1      1    +			 656     677    	24		18		100		32									3
my (@id,@len);
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$id[$tmp[0]]=$tmp[0];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
my (@reads_norm,@rpm);
while (<INPUT>) {
#chrom0 Start1 End2				name3			reads4 strand5 chrom6	Start7		End8 name9  score10 strand11 overlap12
#1       141012  141032  TATGTGTGTATCTTTCGGAAC   2       -       1       140907  141092  1       .       -       20
	chomp;
	my @tmp=split/\t/;
	$reads_norm[$tmp[9]]+=$tmp[4]/$loci{$tmp[3]};
}
close INPUT;
open OUTPUT,">$opts{o}";
for  (1..$#id) {
	if (defined($reads_norm[$_])) {
		$rpm[$_]=$reads_norm[$_]*1000000/$total;
		printf OUTPUT "$_\t%.2f\n",$rpm[$_];
	}else{
		print OUTPUT "$_\t0\n";
	}
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
Usage:    intersect_bed_rpm.pl -d gff_cluster -i intersect_gff_result -t loci_database -o result
Function: use intersect result get assigned region RPM
Command:  -i inersect_bed result
          -o cluster_rpm
          -d cluster
		  -t loci database
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-09-13
Notes:    
\n/
    )
}