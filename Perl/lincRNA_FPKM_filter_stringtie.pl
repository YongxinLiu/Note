#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:e:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my %expression; #expression in hash
#tracking_id     gene_id AA1_FPKM        AA2_FPKM
#TCONS_00000001  XLOC_000001     0       0
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$expression{$tmp[0]}=$tmp[1];
	for $i (2..$#tmp) {
		$expression{$tmp[0]}=$tmp[$i] if $expression{$tmp[0]}<$tmp[$i];
	}
}
close DATABASE;

open EXON,"<$opts{e}";
#TCONS_00000049	1
#TCONS_00000130	1
my %exon; #exon in hash
while (<EXON>) {
	chomp;
	my @tmp=split/\t/;
	$exon{$tmp[0]}=$tmp[1];
}
close EXON;



###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#Bed0		 1	2		3			4			5		
#C161913094	93	167	TCONS_00000148	XLOC_000148	-
#C162714414	103	176	TCONS_00000196	XLOC_000196	-
open OUTPUT,">$opts{o}";
#Bed	6exon_number 7expression_max
#C162854572	2	261	TCONS_00000207	XLOC_000207	-	2	512.707
#C162856756	0	270	TCONS_00000208	XLOC_000208	-	2	27872.9
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if ($exon{$tmp[3]}==1 && $expression{$tmp[3]}>2) {
		print OUTPUT "$_\t$exon{$tmp[3]}\t$expression{$tmp[3]}\n";
	}elsif ($exon{$tmp[3]}>1 && $expression{$tmp[3]}>0.5){
		print OUTPUT "$_\t$exon{$tmp[3]}\t$expression{$tmp[3]}\n";
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
Usage:    lincRNA_FPKM_filter_stringtie.pl -i merged.strand.novel3.bed(过滤后的新转录本) -o output_file(按表达过滤的结果) -d database(表达数据) -e exon_number(外显子数量) -h header num
Function: Select FPKM>=0.5 for multiply transcript(isoform), FPKM>=2 for single exon trasncript in at least one samples
Command:  -i merged.strand.novel3.bed(过滤后的新转录本) (Must)
          -o output_file(按表达过滤的结果) (Must)
          -d database(表达数据)
		  -e exon_number(外显子数量)
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-05-21
Notes:    
\n/
    )
}