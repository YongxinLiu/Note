#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i});
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{t}=1 unless defined($opts{t});

open DATABASE,"<$opts{i}";
my @list;
my @cor;
while (<DATABASE>) {
	chomp;
	push @list,$_;
}
close DATABASE;
open R,">$opts{o}";
foreach  (@list) {
	print R "$_\t";
}
print R "\n";
for $i(0..($#list-1)) {
	print R "$list[$i]";
	print R "\t"x($i+1);
	for $j(($i+1)..$#list) {
		$opts{i}=$list[$i];
		$opts{j}=$list[$j];
###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}$opts{j}";
#database in hash
my ($total1,$total2)=(0,0);
my (%sample1,%sample2,%all);
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$sample2{$tmp[0]}=$tmp[1];
	$all{$tmp[0]}+=$tmp[1];
	$total2+=$tmp[1];
}
close DATABASE;

open DATABASE,"<$opts{d}$opts{i}";
#database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$sample1{$tmp[0]}=$tmp[1];
	$all{$tmp[0]}+=$tmp[1];
	$total1+=$tmp[1];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">data";
foreach (keys %all) {
	next if $all{$_}<$opts{t};
	$sample1{$_}=0 unless defined($sample1{$_});
	$sample2{$_}=0 unless defined($sample2{$_});
	$sample1{$_}=$sample1{$_}/$total1*1000000;#read to RPM
	$sample2{$_}=$sample2{$_}/$total2*1000000;
	next if ($sample1{$_}<$opts{t} && $sample2{$_}<$opts{t});
	printf OUTPUT "%.2f\t%.2f\n",$sample1{$_},$sample2{$_};
}
close OUTPUT;
	$R2=`Rscript ~/bin/cor.r`;
	$tmp=(split('\s+',$R2))[1];
	print R "$tmp\t";
	$cor[$i][$j]=$tmp;
	}
	print R "\n";
}
`rm data`;

# sum each samples cor, then sort find the biggest
my @sum;
for $i(0..$#list) {
	for $j(($i+1)..$#list) {
		$sum[$i]+=$cor[$i][$j];
	}
	for $m(0..($i-1)) {
		$sum[$i]+=$cor[$m][$i];
	}
	print R "$list[$i]\t$sum[$i]\n";
}
my ($max,$id)=(0,0);
for  $m(0..$#list) {
	if ($sum[$m]>$max){
		$max=$sum[$m];
		$id=$m;
	}
}
print R "\n$list[$id]\t$sum[$id]\tBest\n\n";
#Find related
for $j(($id+1)..$#list) {
	print R "$list[$j]\t$cor[$id][$j]\n";
}
for $n(0..($id-1)) {
	print R "$list[$n]\t$cor[$n][$id]\n";
}
close R;


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
Usage:    compare_replicate_sRNA3.pl -i sample_list -d sample_directory -o output_file -t threshold of RPM
Function: Calculate multiple samples double-double corelation coeffection. Find the best relation sample, and show the samples relation with it.
		  Procedure: Filter RPM > threshold sRNA; output in temp file data; call R calculate R2; save all the R2 in matrix;
Command:  -i samples list (Must)
		  -o output file, use test
		  -t threshold RPM, default >= 1
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-05-03
Notes:    Change reads to RPM as threshold
		  Batch calculate sRNA correlation coefficient base on the input list.
\n/
    )
}
