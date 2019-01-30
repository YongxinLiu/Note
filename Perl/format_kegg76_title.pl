#!/usr/bin/perl -w
# 加载时间管理，参数管理，文件名和路径处理的基础包，无须安装
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Scripts usage and about.
# 程序的帮助文档，良好的描述是程序重用和共享的基础，也是程序升级和更新的前提
###############################################################################
sub usage {
    die(
        qq!
Usage:    format_kegg76_title.pl -i kegg_all_clean.fa -o kegg_gene_ko_description.txt
Function: Get KgeneID, KO and Kdescription from KEGG database kegg_all_clean.fa 提取基因ID、KO编号和描述
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, liuyongxin_bio\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/9/28
Notes:    
\n!
    )
}

###############################################################################
#命令行参数据的定义和获取，记录程序初始时间，设置参数默认值
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
# 调置参数的初始值，可以添加更多参数的默认值
$opts{h}=1 unless defined($opts{h});

###############################################################################
#读入的数据或注释文件，用于与输入文件比较或注释(可选)，提供三种方式
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
# 1. 散列结构数据库，要求数据文件有唯一ID并且无顺序要求
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}
# 2. 数组结构数据库，无唯一ID，但有顺序要求
#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;
# 3. 批量数据文件，读取一批有相似结构的文件
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
# 正文部分，读取输入文件，列出输入和输入文件的三行作为示例，方便编程处理数据
open INPUT,"<$opts{i}";
#>hsa:163071  ZNF114; zinc finger protein 114; K09228 KRAB domain-containing zinc finger protein
#MLENSRNLAFIDWATPCKTKDATPQPDILPKRTFPEANRVCLTSISSQHSTLREDWRCPK
#TEEPHRQGVNNVKPPAVAPEKDESPVSICEDHEMRNHSKPTCRLVPSQGDSIRQCILTRD
#SSIFKYNPVLNDSQKTHENNEDDGVLGWNIQWVPCGRKTELKSSTWTGSQNTVHHIRDEI
#DTGANRHQRNPFGKAFREDGSLRAHNTHGREKMYDFTQCENTSRNNSIHAMQMQLYTAET
#NKKDCQTGATSANAPNSGSHKSHCTGEKTHKCPECGRAFFYQSFLMRHMKIHTGEKPYEC
#GKCGKAFRYSLHLNKHLRKHVVQKKPYECEECGKVIRESSKYTHIRSHTGEKPYKCKTCG
#KDFAKSSGLKKHLKTHKDEKPCE
#>hsa:2952  GSTT1; glutathione S-transferase theta 1 (EC:2.5.1.18); K00799 glutathione S-transferase [EC:2.5.1.18]
#MGLELYLDLLSQPCRAVYIFAKKNDIPFELRIVDLIKGQHLSDAFAQVNPLKKVPALKDG
open OUTPUT,">$opts{o}";
#hsa:163071	K09228	KRAB domain-containing zinc finger protein
#hsa:2952	K00799	glutathione S-transferase [EC:2.5.1.18]

my %count;
# h参数用于去除有文件头的行
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
	# 可选，输出文件也保留文件头
	#print OUTPUT $tmp;
}
# 输入和输入处理部分，常用按行读取处理并输入，默认按tab分割数据
print OUTPUT "KgeneID\tKO\tKdescription\n";
while (<INPUT>) {
	if (/^>/) {
		chomp;
		$_=~/^>([^ ]+).*(K\d{5}) (.*$)/;
		print OUTPUT "$1\t$2\t$3\n";
	}else{
		next;
	}
}
close INPUT;
close OUTPUT;

###############################################################################
#Record the program running time!
# 输出程序运行时间
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

