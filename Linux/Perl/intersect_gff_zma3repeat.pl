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
#Main text.
###############################################################################
open ANNO,"<$opts{d}";
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#Chr1    omicslab        sRNA_cluster    44570   44956   26.76   -       .       ID=SRC5;length_type=24-nt
#Chr9	source10 type11	start12	end13	score14	strand15 CDS_phased16	detail17(id,name,biotype)	overlap18
#3	.	repeat_region	218848417	218849555	.	+	.	Name=RLC_opie_AC197691_5727;class=LTR/opie;repeat_consensus=N;type=LTRs	2
open OUTPUT,">$opts{o}";
#ID|chr|start|end|strand|overlap|source|type|description;
#AT1G01010|Chr1|3631|5899|+|2605|TAIR10|protein_coding_gene|NAC domain containing protein 1;
my (%anno,%overlap);
while (<ANNO>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+)/;
	my $src_id=$1;
	$overlap{$src_id}=0 unless defined($overlap{$src_id});
	if ($tmp[18]>$overlap{$src_id}) {
		$overlap{$src_id}=$tmp[18];
		$tmp[17]=~/Name=([^;]+);class=([^;]+)/;
		$anno{$src_id}="$2\|$1\|$tmp[9]\|$tmp[12]\|$tmp[13]\|$tmp[15]\|$tmp[18]";
	}
}
close ANNO;

open INPUT,"<$opts{i}";
while (<INPUT>) {
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#1       omicslab        sRNA_cluster    3758    3861    6.57    -       .       ID=SRC1;length_type=24-nt
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+)/;
	$src_id=$1;
	if (!defined($anno{$src_id})) {
		print OUTPUT "$_\tintergenic\n";
	}else{
		print OUTPUT "$_\t$anno{$src_id}\n";
	}
}
close ANNO;
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
Usage:    intersect_gff_zma3repeat.pl -i SRCs -d sect_SRC_anno -o SRC_anno
Function: Show the intersect file as a SRCs and annotation file.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d intersect_SRC_anno
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-07-04
Notes:    Improve recognize all the SRCs ID
\n/
    )
}
