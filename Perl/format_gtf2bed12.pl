#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

##database in array
#my (@tmp1,@tmp2);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;


#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#0chr	1				2type	3start	4end	5		6strand	7			8ID
#5       protein_coding  exon    1579    3920    .       -       .        gene_id "GRMZM2G356204"; transcript_id "GRMZM2G356204_T01"; exon_number "1"; seqedit "false";
#5       protein_coding  CDS     1684    3903    .       -       0        gene_id "GRMZM2G356204"; transcript_id "GRMZM2G356204_T01"; exon_number "1"; protein_id "GRMZM2G356204_P01";
#5       protein_coding  start_codon     3901    3903    .       -       0        gene_id "GRMZM2G356204"; transcript_id "GRMZM2G356204_T01"; exon_number "1";
#5       protein_coding  stop_codon      1681    1683    .       -       0        gene_id "GRMZM2G356204"; transcript_id "GRMZM2G356204_T01"; exon_number "1";
open OUTPUT,">$opts{o}";
#0chr		1min			2max		3transcript_id	4		5strand	6start_min		7stop_min		8		9exon	10exon length											11exon start
#chr1    48998526        50489626        NM_032785       0       -       48999844        50489468        0       14      1439,27,97,163,153,112,115,90,40,217,95,125,123,192,    0,2035,6787,54149,57978,101638,120482,130297,334336,512729,712915,1164458,1318541,1490908,
#chr1    16767166        16786584        NM_018090       0       +       16767256        16785385        0       8       182,101,105,82,109,178,76,1248, 0,2960,7198,7388,8421,11166,15146,18170,
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}

# Initial data, use first line 
my ($chr,$min,$max,$s_min,$s_max,$exon_length,$exon_start,$strand,$exon_num);
my @exon_start=();
my $pointer=0;
$first=<INPUT>;
chomp($first);
my @tmp=split/\t/,$first;
$tmp[8]=~/transcript\_id\s\"([^\"]+)/;
#print $1,"\n";
$chr=$tmp[0];
$min=$tmp[3]-1;
$max=$tmp[4]-1;
$id=$1;
$strand=$tmp[6];
$exon_num=1;
$exon_length=$max-$min+1;
$exon_length.=",";
if ($strand eq '+') {
	$exon_start="0,";
}else{
	push (@exon_start,$tmp[3]);
	$exon_start="";
}
#print $exon_length,"\n",$exon_start,"\n";

# Cycle each line
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[2]=~/(exon|codon|CDS)/;
	next if ($tmp[2]=~/CDS/ && $pointer==1); #有的正向转录本没有填写起始codon，要在CDS中读取
	$tmp[8]=~/transcript\_id\s\"([^\"]+)/;
#	print $1,"\n";
	$tmpid=$1;

	# if in different gene, output and reinitial
	if ($tmpid ne $id) {
		if ($strand eq '-') {
			@exon_start=reverse(@exon_start);
			for $i(0..$#exon_start) {
				$start=$exon_start[$i]-$min-1;
				$exon_start.="$start\,";
#				print $exon_start[$i],"\t",$min,"\t",$start,"\t",$exon_start,"\n"; #采用stdin暂停调试程序输出结果
#				$count++;
#				<stdin> if $count>30;
			}
		}
#		print "$chr\t$min\t$max\t$id\t0\t$strand\t$s_min\t$s_max\t0\t$exon_num\t$exon_length\t$exon_start\n";
		$s_min=$min unless defined($s_min);
		$s_max=$max unless defined($s_max);
		print OUTPUT "$chr\t$min\t$max\t$id\t0\t$strand\t$s_min\t$s_max\t0\t$exon_num\t$exon_length\t$exon_start\n";
		undef($s_min);
		undef($s_max);
		$id=$tmpid;
		$chr=$tmp[0];
		$min=$tmp[3]-1;
		$max=$tmp[4]-1;
		$id=$tmpid;
		$strand=$tmp[6];
		$exon_num=1;
		$exon_length=$max-$min+1;
		$exon_length.=",";
		$pointer=0;
		@exon_start=();
		if ($strand eq '+') {
			$exon_start="0,";
		}else{
			push (@exon_start,$tmp[3]);
		$exon_start="";
		}
	# In same gene
	}else{
		if ($tmp[6] eq '+') {
			if ($tmp[2]=~/CDS/) {
				$s_min=$tmp[3]-1; #从CDS从读取start_codon min
				$pointer=1;
			}
			if ($tmp[2]=~/start\_codon/) {
				$s_min=$tmp[3]-1; #标注编码区起始为start_codon start
			}
			if ($tmp[2]=~/stop\_codon/) {
				$s_max=$tmp[3]-1; #标注编码区起始为end_codon start,因为end_codon不编码
			}
			if ($tmp[2]=~/exon/) {
				$max=$tmp[4]-1;
				$exon_num++;
				$len=$tmp[4]-$tmp[3]+1;
				$exon_length.="$len\,";
				$start=$tmp[3]-$min-1;
				$exon_start.="$start\,";
			}
		}
		if ($tmp[6] eq '-') { # exon顺序仍然按基因组从左到右，距离也是到基因组左起点
			if ($tmp[2]=~/start\_codon/) {
				$s_max=$tmp[4]-1;
			}
			if ($tmp[2]=~/stop\_codon/) {
				$s_min=$tmp[4]-1;
			}
			if ($tmp[2]=~/exon/) {
				$min=$tmp[3]-1;
				$exon_num++;
				$len=$tmp[4]-$tmp[3]+1;
				$exon_length="$len\,".$exon_length; # exon length from chr left to right
#				$exon_length.="$len\,"; # exon length from chr right to left
#				$start=$max-$tmp[4];
#				$exon_start.="$start\,"; #现在负链是到右侧的距离，参考例子应该为到左侧的距离
				push (@exon_start,$tmp[3]);#暂时没有读到左端最小值，保存数组，输出前计算
			}
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
Usage:    format_gtf2bed12.pl -i gtf -o bed12 -h header num
Function: format gtf to bed12 (for IGV, RSeQC); multiple exon into column 10,11; start and end codon into column 6, 7;
Command:  -i inpute gtf file name (Must)
          -o output bed12 file name (Must)
          -d database file name
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-07-28
Notes:    gtf start chr position is 1, but bed start is 0
		  目前用IGV比较bed与gtf基本一致，只有负链的start_codon起始少2bp，而文件数据位置又没有错？目前结果不影响intron, exon结果。
		  对无复杂结构的单exon基因前两行会报末定义变量，
		  对于无标注start and end codon的基因，会标错误的位置，
\n/
    )
}