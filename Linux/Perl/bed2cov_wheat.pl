#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:p:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{s}=1000000 unless defined($opts{s});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{p}";
#chr	length(unit_in_bp)
#1A      123306033
#2A      251497053
my %db; #database in hash
my $full_len;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
#	$db{$tmp[0]}{start}=$tmp[1];
#	$db{$tmp[0]}{end}=$tmp[2];
	$db{$tmp[0]}{len}=int($tmp[1]/$opts{s});
   $full_len+=$tmp[1];
}

close DATABASE;


###############################################################################
#Main text.
###############################################################################
my(@hang,$site,%num,$Bd,$n,$start,$end,$m);
%num=();
 open F1,"$opts{i}";   ### bed-format
 while(<F1>){
	chomp;
   @hang=split;
   $Bd=$hang[0];
   $site=int($hang[1]/$opts{s});
   $num{$Bd}{$site}++;
 } 
close F1;
$average=1000000/$full_len*$opts{s}/2;
open OUT,">$opts{o}";
foreach $m(sort keys %num) { 
  for($n=0;$n<$db{$m}{len};$n++) # only < means remove incomplete bin 
  {$start=$n*$opts{s};
   $end=($n+1)*$opts{s};
#	next if $n<$db{$m}{start} || $n>$db{$m}{end};
	if(exists $num{$m}{$n}){
		$num{$m}{$n}=$num{$m}{$n}/$opts{d}*1000000;
		$fold=$num{$m}{$n}/$average;
	    printf OUT "$m\t$start\t$end\t%.2f\t%.2f\n",$num{$m}{$n},$fold;
	}else{
		print OUT "$m\t$start\t$end\t0\t0\n";
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
Usage:    bed2cov_wheat.pl -i mapped_bed(frombam) -o normalized_rpm_1k -d total_reads -p centromere_position
Function: Calculate mapping reads to 1M RPM, position unit bp
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d total reads
		  -p position include chr start end chr_length
          -h header line number, default 0
          -s bin size, default 1M
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-10-15
Notes:    add copies, remove incomplete bin
\n/
    )
}