#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{s}="TAIR" if (!defined($opts{s}) || $opts{s} eq "ath");
$opts{s}="ensembl" if ($opts{s} eq "zma");


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my ($pointer,$id,$length_type,$distance);
while (<INPUT>) {
#0chr 1	2	3start	4end	5RPM	6strand	7	8SRCID,length	9overlap	10leaft	11right
#1       omicslab        SRC     4261    4364    1.82    +       .       ID=zma-SRC1;length_type=24-nt   intergenic      ;zma-mtec4|1|3619|3931|-|329|MTEC|RLX&Unknown&loukuv|RLX_loukuv_AC197842_5799   GRMZM2G059865|1|4854|9652|-|489|ensembl|protein_coding|AT1G75950&S phase kinase-associated protein 1;zma-mtec5|1|10872|15310|-|6507|MTEC|DTC&CACTA&NA|ZM_CACTA_42
#1       omicslab        SRC     111878  112015  2.00    +       .       ID=zma-SRC4;length_type=24-nt   zma-mtec53|1|111814|111994|-|117|MTEC|RLX&Unknown&small|RLX_small_AC195943_142;zma-mtec54|1|111872|113431|-|138|MTEC|RLX&Unknown&small|RLX_small_AC217574_13522;        GRMZM2G093344|1|109519|111769|-|108|ensembl|protein_coding|LOC_Os01g09090&expressed protein;zma-mtec52|1|111146|111473|-|404|MTEC|RLC&Copia&guvi|RLC_guvi_AC185473_1128 GRMZM2G093399|1|136307|138929|+|24291|ensembl|protein_coding|NA;zma-mtec55|1|112152|114976|-|136|MTEC|RLX&Unknown&ebel|RLX_ebel_AC188777_2128
#1       omicslab        SRC     234920  234999  8.67    +       .       ID=zma-SRC10;length_type=other  GRMZM2G104572|1|234779|235890|+|80|ensembl|protein_coding|NA;   AC177838.2_FG015|1|161143|161925|+|72994|ensembl|protein_coding|NA;zma-mtec146|1|234300|234434|-|485|MTEC|RLC&Copia&dijap|RLC_dijap_AC211466_11198      GRMZM5G822187|1|243464|243742|+|8464|ensembl|protein_coding|NA;zma-mtec147|1|246009|246104|-|11009|MTEC|RLG&Gypsy&dagaf|RLG_dagaf_AC195302_4533
	chomp;
	$pointer=0;#pointer, if one SRCs overlapped genes, not calculate distance to nearest gene.
	$distance=0;
	$geneID="";
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+);length_type=([\w\-]+)/;
#	print "$1\t$2\n";
	$id=$1;
	$length_type=$2;
#	print "$tmp[9]\n";
	my @overlap=split(';',$tmp[9]);
	foreach  my $tmp2(@overlap) {
		if ($tmp2=~/$opts{s}/) {
			#print "$tmp2\n";
			my @tmp3=split('\|',$tmp2);
			#0id	1chr	2start	3end	4strand	5overlap	6	7type	8description
			#GRMZM2G104572|1|234779|235890|+|80|ensembl|protein_coding|NA
			$distance=(($tmp[3]+$tmp[4])/2-$tmp3[2]+1)/($tmp3[3]-$tmp3[2]+1);
			#print "$distance\n";
			#May some SRCs partial overlap genes, if overlap > 1, set to 0.99, if overlap <0, set to 0.01, wrong
			#May some SRCs partial overlap genes, if overlap > 1,or overlap <0, calculate distance to gene
			if ($distance>0 && $distance<1) {
				$distance=1-$distance if $tmp3[4] eq '-';
			}
			if ($distance>=1 || $distance <=0) {
				$distance=($tmp[3]+$tmp[4])/2-$tmp3[3]+1;
				$distance*=(-1) if $tmp3[4] eq '-';
			} 
			$geneID=$tmp3[0];
			#print "$distance\n";
			$pointer=1;
			last;
		}
	}
	
	if ($pointer==0){
		my @left=split(';',$tmp[10]);
		#print "$tmp[10]\n";
		foreach  my $tmp2(@left) {
			#GRMZM2G059865|1|4854|9652|-|489|ensembl|protein_coding|AT1G75950&S phase kinase-associated protein 1
			if ($tmp2=~/$opts{s}/) {
				my @tmp3=split('\|',$tmp2);
				$distance=$tmp3[5];
				$distance=$distance*(-1) if $tmp3[4] eq '-';
				$geneID=$tmp3[0];
			}
		}
		#print "$distance\n";
		my @right=split(';',$tmp[11]);
		#print "$tmp[11]\n";
		foreach  my $tmp2(@right) {
			#GRMZM2G059865|1|4854|9652|-|489|ensembl|protein_coding|AT1G75950&S phase kinase-associated protein 1
			if ($tmp2=~/$opts{s}/) {
				#print "$tmp2\n";
				my @tmp3=split('\|',$tmp2);
				#print "$tmp3[5]\n";
				if ($tmp3[5]<abs($distance) || $distance==0) {
					$distance=$tmp3[5];
					$distance=$distance*(-1) if $tmp3[4] eq '+';
					$geneID=$tmp3[0];
				}
			}
		}
	}

	next unless defined($distance);
	if ($geneID=~/\.\d$/) {
		$geneID=(split('\.',$geneID))[0];
	}
	if ($geneID=~/\-\d$/) {
		$geneID=(split('\-',$geneID))[0];
	}
	print OUTPUT "$id\t$length_type\t$distance\t$geneID\n";

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
Usage:    annotation_gene_distribution.pl -i SRCs_overlap_flank -o relation_with_gene -d database	-h header num
Function: Based on SRCs overlap and flank gene annotation, get the SRCs distribution around genes.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
		  -s species select, default ath, other include zma
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.2
Update:   2014-10-11
Notes:    v1.1 Debug overlap region value >1 or < 0, but have peak in gene start and end
		  v1.2 overlap region value >1 or < 0 all classify into up or down
\n/
    )
}