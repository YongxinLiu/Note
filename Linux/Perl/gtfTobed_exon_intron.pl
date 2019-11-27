#!/usr/bin/perl
use warnings;
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################

open INPUT,"<$opts{i}";
#0Chr		1gene		2type	3start	4end	5	6strand	7	8detail
#scaffold262259	protein_coding	exon	21521	21569	.	-	.	 gene_id "TRIUR3_26600"; transcript_id "TRIUR3_26600-T1"; exon_number "1"; seqedit "false";
#scaffold262259	protein_coding	CDS	21521	21569	.	-	0	 gene_id "TRIUR3_26600"; transcript_id "TRIUR3_26600-T1"; exon_number "1"; protein_id "TRIUR3_26600-P1";
#scaffold262259	protein_coding	start_codon	21567	21569	.	-	0	 gene_id "TRIUR3_26600"; transcript_id "TRIUR3_26600-T1"; exon_number "1";
open OUTPUT,">$opts{o}";
#chr0			start1  end2              ID3       value4    strand5
#scaffold262259    8       31      TRIUR3_26600        exon/intron       -
$pointer=""; # set pointer, when id same, write intron
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[2] eq "exon";
#	print $tmp[8],"\n";
	$tmp[8]=~/gene_id\s\"([\w\_]+)\"/;
	my $gene=$1;
#	print "$1\n";
	if ($gene eq $pointer) {
		if ($tmp[6] eq '+') {
			$end=$tmp[3]-1;
			print OUTPUT "$tmp[0]\t$start\t$end\t$gene\tintron\t$tmp[6]\n";
			$start=$tmp[4]+1;
		}else{
			$start=$tmp[4]+1;
			print OUTPUT "$tmp[0]\t$start\t$end\t$gene\tintron\t$tmp[6]\n";
			$end=$tmp[3]-1;
		}
		print OUTPUT "$tmp[0]\t$tmp[3]\t$tmp[4]\t$gene\texon\t$tmp[6]\n";
	}else{
		print OUTPUT "$tmp[0]\t$tmp[3]\t$tmp[4]\t$gene\texon\t$tmp[6]\n";
		$pointer=$gene;
		if ($tmp[6] eq '+') {
			$start=$tmp[4]+1;
		}else{
			$end=$tmp[3]-1;
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
Usage:    gtfTobed_exon_intron.pl -i gtf file (exon, cds, start-codon ...) -o gtf file (exon and intron)
Function: get standard gtf file, and output only contain exon and intron lines
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-07-14
Notes:    
\n/
    )
}