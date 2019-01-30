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
#0Chr		1gene		2type	3start	4end	5	6strand	7	8detail
#scaffold125401  ensembl miRNA_gene      59      183     .       -       .       ID=gene:EPlTURG00000000616;biotype=miRNA;description=microRNA MIR444 [Source:RFAM;Acc:RF00920];external_name=MIR444;logic_name=ncrna_eg
#scaffold125401  ensembl miRNA   59      183     .       -       .       ID=transcript:EPlTURT00000000616;Parent=gene:EPlTURG00000000616;biotype=miRNA;external_name=MIR444-201;logic_name=ncrna_eg
#scaffold125401  ensembl miRNA_gene      59      183     .       +       .       ID=gene:EPlTURG00000000747;biotype=miRNA;description=microRNA MIR444 [Source:RFAM;Acc
#3       rep14   repeat  218914213       218915768       7076    +       .       ID=ath-rep642857;Note=LTR/Gypsy;Annotation=HUCK1-LTR_ZM
#3       rep14   repeat  218921009       218921154       454     +       .       ID=ath-rep642864;Note=LTR/Gypsy;Annotation=Gypsy27-ZM_LTR
#3       rep14   repeat  218921155       218922776       8594    +       .       ID=ath-rep642865;Note=LTR/Gypsy;Annotation=HUCK1-LTR_ZM
open OUTPUT,">$opts{o}";
#chr0			start1  end2              ID3       value4    strand5
#scaffold262259    8       31      TRIUR3_26600        gene/upstream/downstream       -
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[2]=~/repeat/;
	$tmp[8]=~/ID=([^;]+)/;
	my $gene=$1;
	if ($tmp[6] eq '+') {
		$start=$tmp[3]-$opts{s};
		$start=1 if $start<1;
		print OUTPUT "$tmp[0]\t$start\t",$tmp[3]-1,"\t$gene\tupstream\t$tmp[6]\n" if $tmp[3]-1>=$start;
		print OUTPUT "$tmp[0]\t$tmp[3]\t$tmp[4]\t$gene\tgene\t$tmp[6]\n";
		print OUTPUT "$tmp[0]\t",$tmp[4]+1,"\t",$tmp[4]+$opts{s},"\t$gene\tdownstream\t$tmp[6]\n";
	}else{
		$start=$tmp[3]-$opts{s};
		$start=1 if $start<1;
		print OUTPUT "$tmp[0]\t",$tmp[4]+1,"\t",$tmp[4]+$opts{s},"\t$gene\tupstream\t$tmp[6]\n";
		print OUTPUT "$tmp[0]\t$tmp[3]\t$tmp[4]\t$gene\trepeat\t$tmp[6]\n";
		print OUTPUT "$tmp[0]\t$start\t",$tmp[3]-1,"\t$gene\tdownstream\t$tmp[6]\n" if $tmp[3]-1>=$start;;
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
Usage:    gffTobed_repeat_flank.pl -i gff3 file (repeat, gene) -o bed file (exon and intron)
Function: input standard gff3 file, and output only gene and flank in bed
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
          -s distance, default = 2000
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-07-14
Notes:    
\n/
    )
}