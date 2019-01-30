#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq!
Usage:    format_list2plate.pl -i L1_split.count -o output_dir -d database	-h header num
Function: format each pole into 96 wells plate and 100 wells gels
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/18
Notes:    
\n!
    )
}

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
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
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
my %plate; # 保存所有库板ID
my %well; # 保存所有孔ID及数据
open INPUT,"<$opts{i}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[0]=~/(L\d+P\d+)/;
	$plate{$1}++;
	$well{$tmp[0]}=$tmp[1];
}
close INPUT;

`mkdir -p $opts{o}`;
@j=qw(A B C D E F G H);
foreach $plate (keys %plate) {
open OUTPUT,">$opts{o}$plate.plate";
# 循环板中的96孔
print OUTPUT "Plate\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\n";
foreach $j (@j) {
	print OUTPUT "$j";
	for $i(1..12) {
		$id="$plate$j$i";
		print OUTPUT "\t$well{$id}";
	}
	print OUTPUT "\n";
}
close OUTPUT;

open OUTPUT,">$opts{o}$plate.gel";
# 循环板中的96孔
print OUTPUT "Gel\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\t23\t24\n";
foreach $j (0..3) {
	print OUTPUT $j+1;
	$j=$j*2;
	$k=$j+1;
	for $i(1..12) {
		$id="$plate$j[$j]$i";
#		print $j,"\t",$id,"\n";
		print OUTPUT "\t$well{$id}";
		$id="$plate$j[$k]$i";
		print OUTPUT "\t$well{$id}";
	}
	print OUTPUT "\n";
}
close OUTPUT;
}
###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

