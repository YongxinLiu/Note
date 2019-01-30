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
my(@hang,$site,%num,$Bd,%site,$n,$start,$end,$m);
%num=();

 open F1,"$opts{i}";   ### bed-format
 while(<F1>){
	chomp;
   @hang=split;
   $Bd=$hang[0];
   $site=int($hang[1]/10000);
   $num{$Bd}{$site}++;
 } 
close F1;

open OUT,">$opts{o}";
$site{"Bd1"}=7511;
$site{"Bd2"}=5938;
$site{"Bd3"}=5993;
$site{"Bd4"}=4877;
$site{"Bd5"}=2856;

for($m=1;$m<=5;$m++)
 {$Bd="Bd".$m;
 for($n=0;$n<=$site{$Bd};$n++)
  {$start=$n/100;
   $end=($n+1)/100;
	if(exists $num{$Bd}{$n}){
		$num{$Bd}{$n}=$num{$Bd}{$n}/$opts{d}*1000000;
	   printf OUT "$Bd\t$start\t$end\t%.2f\n",$num{$Bd}{$n};
	}else{
		print OUT "$Bd\t$start\t$end\t0\n";
	}
  }
 }
close OUT;

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
Usage:    bed2_10k-cov_bdi.pl -i mapped_bed(frombam) -o normalized_rpm_10k -d total_reads -h header num
Function: Calculate mapping reads to 10kb RPM, position unit Mb
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-01-26
Notes:    
\n/
    )
}