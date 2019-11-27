#!/usr/bin/perl

use strict;
use Getopt::Std;
use vars qw($opt_i $opt_x $opt_y $opt_l $opt_f $opt_o $opt_h);
getopts('i:x:y:l:f:o:h:');
my $fq_file        = $opt_i;
my $adapter5       = $opt_x ? $opt_x : "GTTCAGAGTTCTACAGTCCGACGATC";
my $adapter3       = $opt_y ? $opt_y : "TCGTATGCCGTCTTCTGCTTG";
my $min_len        = $opt_l ? $opt_l : 18;
my $format         = $opt_f ? $opt_f : 1;
my $out_file       = $opt_o;
my $help           = $opt_h ? 1 : 0;

my $usage = << "USAGE";

Description: Perl script used to filter low quality short reads, remove polyA and trim 3' 5' adapter

Usage: perl format_fastq.pl [options] > clean.txt
Options:
  -i <file>  Short reads file in fastq format 
  -x <str>   5\' adaptor sequence, default="GTTCAGAGTTCTACAGTCCGACGATC"
  -y <str>   3\' adaptor sequence, default="TCGTATGCCGTCTTCTGCTTG"
  -l <int>   The minmal length of the reads. default=18
  -f <int>   Fastq file format: 1=Sanger format; 2=Solexa/Illumina 1.0 format; 3=Illumina 1.3+ format; default=1
  -o <str>   Output file
  -h         Help
Examples: perl format_fastq.pl -i sample.fq -o clean.txt
	  or:
          perl format_fastq.pl -i sample.fq -x "GTTCAGAGTTCTACAGTCCGACGATC" -y "TCGTATGCCGTCTTCTGCTTG" -l 18 -f 2 -o clean.txt
USAGE


if ($help) {
	print $usage;
	exit;
}

unless (( -e $fq_file ) and ($out_file) and ( $adapter5 =~ /[A|T|C|G|U]/i ) and ( $adapter3 =~ /[A|T|C|G|U]/i ) and ( $format =~ /^[1|2|3]$/ )) {
	print $usage;
	exit;
}

my %reads_number;

open IN,"$fq_file" or die $!;
open OUT, ">$out_file" or die $!;
open OUT1, ">report.txt" or die $!;

while (<IN>) { #while读取4行循环的fastq,一行fastq ID,再取序列，再扔一行，再取质量
	my $read = <IN>;
	chomp $read;
	<IN>;
	my $quality_line = <IN>;
	chomp $quality_line;
	if ($read !~ /^[A|T|C|G|U]+$/i) {
		next;
	} else {
		$read =~ tr/atcgu/ATCGU/;
		$adapter3 =~ tr/atcgu/ATCGU/;
		$adapter5 =~ tr/atcgu/ATCGU/;
		my $quality = 0; #对不同格式质量控制，分支很流畅
		if ($format == 2) {
			$quality= &check_qlt($quality_line,64,9);
		}
		elsif ($format == 1) {
			$quality= &check_qlt($quality_line,33,15);
		}
		elsif ($format == 3) {
			$quality= &check_qlt($quality_line,64,10);
		}
		if ($quality == 1) {
			$reads_number{$read}++;
		}
	
	}
}
close IN;

###remove 3' adapter
my %after_mv3;
foreach my $read (keys %reads_number) {	
	my $read_af3 = &mvadapter($read,$adapter3,0);
	
	$after_mv3{$read_af3}+=$reads_number{$read} unless(length($read_af3) < 16);
}
undef(%reads_number);

###remove 5' adapter
#my %after_mv5;
#foreach my $read (keys %after_mv3) {	
#	my $read_af5= &mvadapter($read,$adapter5,1);
#	$after_mv5{$read_af5}+=$after_mv3{$read} unless(length($read_af5) < 16);
#}
#undef(%after_mv3);

###remove polyA
#my %clean;
#foreach my $read (keys %after_mv5) {	
#	my $polya= &mvpolyA($read);
#	unless($polya == 1) {
#		$clean{$read}=$after_mv5{$read};
#	}
#}
#undef(%after_mv5);
my %clean;
foreach my $read (keys %after_mv3) {	
	my $polya= &mvpolyA($read);
	unless($polya == 1) {
		$clean{$read}=$after_mv3{$read};
	}
}
undef(%after_mv3);

###print 
my $uniq_seq = 0;
my $redu_seq = 0;
foreach my $r (sort {$clean{$b} <=> $clean{$a}} keys %clean ) {
	$uniq_seq ++;
	$redu_seq += $clean{$r};
	print OUT "$r\t$clean{$r}\n";
}
print OUT1 "Input file: $fq_file\n";
print OUT1 "Output file: $out_file\n";
print OUT1 "5' adaptor: $adapter5\n";
print OUT1 "3' adaptor: $adapter3\n";
print OUT1 "The minmal length of the reads: $min_len\n";
print OUT1 "This dataset contains $redu_seq redundant sequence reads which can be merged to $uniq_seq non-redundant sequences.\n";

sub check_qlt
{
	my $quality_line = shift;
	my $asc = shift;
	my $tv = shift;
	my $num = 0;
	my $count = 0;
	my @ql = split (//,$quality_line);
	my $wid = $#ql+1;
	foreach my $i (0..$#ql)
	{
		$num = ord($ql[$i])-$asc;
		if($num-$tv <0){$count++;}
	}

		if( $count > 1 ){return 0;}
		else {return 1;}
	
}


sub mvadapter {
	my $read = shift; #shift直接取默认数组元素，从头开始
	my $adapter = shift;
	my $mode = shift;
	my $readback;
	if ($mode == 1) {
		$read = reverse($read);
		$adapter = reverse($adapter);	
	}
	my $bl= length($read);
	my $tl= length($adapter);
	my @bemapped=split(//,$read);
	my @tomap=split(//,$adapter);
	my @record;
	for (my $i =0; $i<$bl;$i++) {#采用接头在序列上滑动比较法，进行统计
		my $match =0;
		my $mismatch = 0;
		for (my $n=0;$n<$tl;$n++) {
			last unless( $bemapped[$i+$n]); #超出序列长度的则终止
			if($bemapped[$i+$n] eq $tomap[$n]) {#滑动比较法
				$match++;
			} else {
				$mismatch++;
				last if ($mismatch >3);#找到3个以内错配的接送就算匹配成功，大于三个的终止
			}	
		}
		my $long= $match+$mismatch;
		my $per = sprintf "%.2f",$mismatch/($match+$mismatch);
		if($mismatch < 4 and $per < 0.3 and $long >4) { #输出4个以上错配，错配率小于30%且匹配长度在4个以内
			push @record ,[$per,$i,$mismatch,$match];
		}
	}
	
	if ($#record == 0) { 
		$readback = substr($read,0,$record[0][1]);
		if ($mode ==1){$readback = reverse $readback;}
	}
	elsif($#record > 0) {
		my @record_sort = sort {$a->[0] <=> $b->[0]} @record;
		$readback = substr($read,0,$record_sort[0][1]);
		if ($mode ==1){$readback = reverse $readback;}
	} else {
		$readback = $mode ? reverse($read) : "null";
	}
	return $readback;
}

sub mvpolyA {
	my $read = shift;
	$read =~ /(A{3,})/i;
	my $rp=$1;	
	if ((length $rp) > 3) {
		return 1;
	} else {
		return 0;
	}
}
