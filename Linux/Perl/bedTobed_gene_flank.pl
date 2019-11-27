#!/usr/bin/perl
use warnings;
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o});
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{s}=2000 unless defined($opts{s});

###############################################################################
#Main text.
###############################################################################

open INPUT,"<$opts{i}";
#ID0			start1	end2	id3	length4	strand5
#scaffold59710	49899	52893	3	2994	-
#scaffold26969	20735	23729	4	2994	+
#0Chr		1gene		2type	3start	4end	5	6strand	7	8detail
#scaffold125401  ensembl miRNA_gene      59      183     .       -       .       ID=gene:EPlTURG00000000616;biotype=miRNA;description=microRNA MIR444 [Source:RFAM;Acc:RF00920];external_name=MIR444;logic_name=ncrna_eg
#scaffold125401  ensembl miRNA   59      183     .       -       .       ID=transcript:EPlTURT00000000616;Parent=gene:EPlTURG00000000616;biotype=miRNA;external_name=MIR444-201;logic_name=ncrna_eg
#scaffold125401  ensembl miRNA_gene      59      183     .       +       .       ID=gene:EPlTURG00000000747;biotype=miRNA;description=microRNA MIR444 [Source:RFAM;Acc
open OUTPUT,">$opts{o}";
#chr0			start1  end2              ID3       value4    strand5
#scaffold262259    8       31      TRIUR3_26600        gene/upstream/downstream       -
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
#	next unless $tmp[2]=~/gene/;
#	$tmp[8]=~/gene:([^;]+)/;
#	my $gene=$1;
	if (!defined($tmp[5])) {
		$tmp[5]='+';
	}
	if (($tmp[5] eq '+')) {
		$start=$tmp[1]-$opts{s};
		$start=1 if $start<1;
		print OUTPUT "$tmp[0]\t$start\t",$tmp[1]-1,"\t$tmp[3]\tupstream\t$tmp[5]\n" if $tmp[1]-1>=$start;
		print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\tgene\t$tmp[5]\n";
		print OUTPUT "$tmp[0]\t",$tmp[2]+1,"\t",$tmp[2]+$opts{s},"\t$tmp[3]\tdownstream\t$tmp[5]\n";
	}else{
		$start=$tmp[1]-$opts{s};
		$start=1 if $start<1;
		print OUTPUT "$tmp[0]\t",$tmp[2]+1,"\t",$tmp[2]+$opts{s},"\t$tmp[3]\tupstream\t$tmp[5]\n";
		print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\tgene\t$tmp[5]\n";
		print OUTPUT "$tmp[0]\t$start\t",$tmp[1]-1,"\t$tmp[3]\tdownstream\t$tmp[5]\n" if $tmp[1]-1>=$start;;
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
Usage:    bedTobed_gene_flank.pl -i bed file (repeat, gene) -o bed file (exon and intron)
Function: input standard gff3 file, and output only gene and flank in bed
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
          -s distance, default = 2000
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-07-28
Notes:    
\n/
    )
}