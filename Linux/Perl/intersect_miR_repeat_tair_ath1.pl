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
#Chr9	source10 type11	start12	end13	score14	strand15 CDS_phased16	detail17(id,name,biotype)
#Chr1    RepeatMasker    repeat  1       105     .       -       .       ID=1;Note=DNA;Annotation=ATREP18        98
#Chr1    TAIR10  gene    44677   44787   .       +       .       ID=AT1G01073;Note=protein_coding_gene;Annotation=protein_coding_gene    111
open OUTPUT,">$opts{o}";
#ID|chr|start|end|strand|overlap|source|type|description;
#AT1G01010|Chr1|3631|5899|+|2605|TAIR10|protein_coding_gene|NAC domain containing protein 1;
my (%anno,%count);
while (<ANNO>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+)/;
	my $src_id=$1;
#	$anno{$src_id}{dis}=0 unless defined($anno{$src_id}{dis});
	if (/Parent=/) {
		$tmp[17]=~/Parent=([^,^;]+)/;
		$id=$1;
		$id=~s/_T01//; #format maize transcript ID to gene ID
		$id=~s/FGT/FG/;
		$count{$id}++;
		$anno{$src_id}.="$id\-$count{$id}\|$tmp[9]\|$tmp[12]\|$tmp[13]\|$tmp[15]\|$tmp[18]\|$tmp[10]\|$tmp[11]\|$id\;";
	}else{
		$tmp[17]=~/ID=([^;]+);Note=([^;]+);Annotation=([^;]*)/;
		$anno{$src_id}.="$1\|$tmp[9]\|$tmp[12]\|$tmp[13]\|$tmp[15]\|$tmp[18]\|$tmp[10]\|$2\|$3\;";

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
Usage:    intersect_miR_repeat_tair_ath1.pl -i SRCs -d sect_SRC_anno -o SRC_anno
Function: Show the intersect file as a SRCs and annotation file.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d intersect_SRC_anno
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-06
Notes:    Improve recognize all the SRCs ID
\n/
    )
}
