#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:l:r:s:m:n:p:k:x:g:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{l}=20 unless defined($opts{l});
$opts{r}=2 unless defined($opts{r});
$opts{s}=100 unless defined($opts{s});
$opts{m}=3 unless defined($opts{m});
$opts{n}=100 unless defined($opts{n});
$opts{p}=1 unless defined($opts{p});
$opts{k}=5 unless defined($opts{k});
$opts{x}="n" unless defined($opts{x});
print "You select parameter:\nsequences_loci<=$opts{l}\nsequences_reads>=$opts{r}\nsequences_distance<=$opts{s}\ncluster_loci>=$opts{m}\ncluster_reads>=$opts{n}\nRPM>$opts{p}\nRPKM>$opts{k}\nOutput select mapping:$opts{x}\n";

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#sequences0			reads1	loci2
#AAAAAAAAAAAAAACTCGGCC   1       1
my ($total,$reads_in_cluster,$length_in_cluster,$reads_per,$length_per);#toal mapped reads,all reads in cluster, all length in cluster, reads percentage, length percentage
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[2];#preserve all the loci number in hash
	$total+=$tmp[1];
}
close DATABASE;
print "\nDatabase total reads:$total\n";

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#id&reads0	chr1	strand2		start3	end4 	sequences5
#21719662&1      1       +       656     677     TGGATATAAACTATTTTTGGCT
open OUTPUTR,">$opts{o}";
#id0	chr1 strand2	start3	end4 	length5	loci6	reads7	reads_nrom8(normalized,reads/loci)	RPM8	RPKM9(reads_norm/KM)	propotion(normalized reads in dominant strand)
#1      1    +			 656     677    	24		18		100		32								0.21	2.13						0.93
my ($id,$chr,$start,$end,$length,$loci,$reads,$reads_norm,$rpm,$rpkm,$propotion,$reads_watson)=(1,0,0,0,0,0,0,0,0,0,0,0);#reads_watson record the watson strand normalized reads, for compute propotion and decide the dominant strand
my ($len,%length);#length hash, record length
#Use file line initial variation
$tmp=<INPUT>;
chomp($tmp);
@tmp=split/\t/,$tmp;
$tmp[0]=~/(\d+)&(\d+)/;
$chr=$tmp[1];
$start=$tmp[3];
$end=$tmp[4];
$loci=1;
$reads=$2;
$reads_norm=($2/$database{$tmp[5]});
$reads_watson=$reads_norm if $tmp[2] eq "+";
$len=length($tmp[5]);
$length{$len}+=$reads_norm;
while (<INPUT>) {
	chomp;
	@tmp=split/\t/;
	$tmp[0]=~/(\d+)&(\d+)/;#get id in $1 and reads number in $2
	next if $2<$opts{r};#filter low than $opts{r} reads sequences
	next if $database{$tmp[5]}>$opts{l};#filter more than $opts{l} hits sequences
#	if ($opts{x} eq "y") {
#		open OUTPUT,">>$opts{o}.psmap";
#		#id0	reads1	chr2	strand3	start4	end5 	loci6	sequences6
#		#21719662	1      1       +       656     677     18	TGGATATAAACTATTTTTGGCT
#		print OUTPUT "$1\t$2\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$database{$tmp[5]}\t$tmp[5]\n";#input select mapping loci
#		close OUTPUT;
#	}
	if ($chr eq $tmp[1] && $tmp[3]-$end>$opts{s}+1) {#In the same chr and distance more than threshold of distance($opts{s})
		$length=$end-$start+1;
		$rpm=$reads_norm*1000000/$total;
		$rpkm=$reads_norm*1000000*1000/$total/$length;
		$propotion=$reads_watson/$reads_norm;
		if ($propotion<0.5) {
			 $strand="-";
			 $propotion=1-$propotion;
		}else{
			 $strand="+";
		}
		if ($loci>=$opts{m} && $reads>=$opts{n} && $rpm>=$opts{p} && $rpkm>=$opts{k}) {
			printf OUTPUTR "$id\t$chr\t$strand\t$start\t$end\t$length\t$loci\t$reads\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",$reads_norm,$rpm,$rpkm,$length{20},$length{21},$length{22},$length{23},$length{24},$propotion;
			$reads_in_cluster+=$reads_norm;
			$length_in_cluster+=$length;
			$id++;
		}
		$start=$tmp[3];
		$end=$tmp[4];
		$loci=1;
		$reads=$2;
		$reads_norm=($2/$database{$tmp[5]});
		$reads_watson=0;
		$reads_watson=$reads_norm if $tmp[2] eq "+";
		for (20..24) {
			$length{$_}=0;
		}
		$len=length($tmp[5]);
		$length{$len}+=($2/$database{$tmp[5]});
	}elsif ($chr eq $tmp[1] && $tmp[3]-$end<=$opts{s}+1){#In the same cluster
		$start=$tmp[3] if $tmp[3]<$start;
		$end=$tmp[4] if $tmp[4]>$end;
		$loci++;
		$reads+=$2;
		$reads_norm+=($2/$database{$tmp[5]});
		$reads_watson+=($2/$database{$tmp[5]}) if $tmp[2] eq "+";
		$len=length($tmp[5]);
		$length{$len}+=($2/$database{$tmp[5]});
	}elsif ($chr ne $tmp[1]){
		$length=$end-$start+1;
		$rpm=$reads_norm*1000000/$total;
		$rpkm=$reads_norm*1000000*1000/$total/$length;
		$propotion=$reads_watson/$reads_norm;
		if ($propotion<0.5) {
			 $strand="-";
			 $propotion=1-$propotion;
		}else{
			 $strand="+";
		}
		if ($loci>=$opts{m} && $reads>=$opts{n} && $rpm>=$opts{p} && $rpkm>=$opts{k}) {
			printf OUTPUTR "$id\t$chr\t$strand\t$start\t$end\t$length\t$loci\t$reads\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",$reads_norm,$rpm,$rpkm,$length{20},$length{21},$length{22},$length{23},$length{24},$propotion;
			$reads_in_cluster+=$reads_norm;
			$length_in_cluster+=$length;
			$id++;
		}
		$chr=$tmp[1];
		$start=$tmp[3];
		$end=$tmp[4];
		$loci=1;
		$reads=$2;
		$reads_watson=0;
		$reads_norm=($2/$database{$tmp[5]});
		$reads_watson=$reads_norm if $tmp[2] eq "+";
		for (20..24) {
			$length{$_}=0;
		}
		$len=length($tmp[5]);
		$length{$len}+=($2/$database{$tmp[5]});
	}
}
$length=$end-$start+1;
$rpm=$reads_norm*1000000/$total;
$rpkm=$reads_norm*1000000*1000/$total/$length;
$propotion=$reads_watson/$reads_norm;
if ($propotion<0.5) {
	 $strand="-";
	 $propotion=1-$propotion;
}else{
	 $strand="+";
}
if ($loci>=$opts{m} && $reads>=$opts{n} && $rpm>=$opts{p} && $rpkm>=$opts{k}) {
	printf OUTPUTR "$id\t$chr\t$strand\t$start\t$end\t$length\t$loci\t$reads\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",$reads_norm,$rpm,$rpkm,$length{20},$length{21},$length{22},$length{23},$length{24},$propotion;
	$reads_in_cluster+=$reads_norm;
	$length_in_cluster+=$length;
}
close INPUT;
close OUTPUTR;
$reads_per=$reads_in_cluster/$total*100;
$length_per=$length_in_cluster/$opts{g}*100;
printf "These clusters incorporate %.2f%% (%.0f) of total mapped reads.\n",$reads_per,$reads_in_cluster;
printf "These clusters coverage is %.2f%% of genome.\n",$length_per;


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
Usage:    define_cluster6.pl -i inpute_file -d loci_file -o output_file
Function: define sRNA cluster in the genome
		  (1) Sequence need to perfect match reference genome; 
		  (2) loci < 21 mapping sequences for exclude highly repeat sequences; 
		  (3) reads > 1 in order to exclude occasion transcription or contaminate; 
		  (4) any sRNA distance <= 200bp as in one cluster, any cluster distance > 200bp as independence cluster. 
Command:  -i input file, Mapping-v5 mapping result
          -o output file, filtered sRNA mapping result in *.vt, defined cluster summary in *.sum
          -d sequences loci result
          -l threshold of sequences loci, default 20, exclude reads have more than 20 loci on genome
          -r threshold of sequences reads, default 2, exclude reads less than 2 sequences
          -s threshold of distance, default 200, any sRNA distance < 200bp as in one cluster, any cluster distance > 200bp as independence cluster.
          -m threshold of cluster loci number, default 3, exclude cluster less than 3 loci
          -n threshold of cluster reads, default 100
          -p threshold of reads per million (RPM, normalized expression), default 1
          -k threshold of reads per kilobase of million (RPKM, normalized expression), default 5
          -x output select mapping result y\/n, default n
          -g genome size, for calculate genome coverage e.g. Ath 121182535 Zma 2100873268
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v6.0
Update:   2013-05-02
Notes:    v2.0 add threshold of loci, reads, RPKM
		  v2.1 add sRNA used and genome coverage
		  v3.0 add each cluster sRNA length 20-24 count
		  v4.0 change RPKM to RPM, not good
		  v5.0 cluster not seperate strand, but add preference strand, and prefer rate
		  v6.0 add RPM before RPKM, add RPM threshold
\n/
    )
}