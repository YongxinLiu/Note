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
$opts{s}=10000 unless defined($opts{s});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{p}";
#chr	start	end	length	(unit in bp)
#1	134124000	135251000	301476924
my %db; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$db{$tmp[0]}{start}=$tmp[1];
	$db{$tmp[0]}{end}=$tmp[2];
	$db{$tmp[0]}{len}=int($tmp[3]/$opts{s});
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

open OUT,">$opts{o}";
for($m=1;$m<=10;$m++)
 {for($n=0;$n<=$db{$m}{len};$n++)
  {$start=$n/100;
   $end=($n+1)/100;
#	next if $n<$db{$m}{start} || $n>$db{$m}{end};
	if(exists $num{$m}{$n}){
		$num{$m}{$n}=$num{$m}{$n}/$opts{d}*1000000;
	   printf OUT "$m\t$start\t$end\t%.2f\n",$num{$m}{$n};
	}else{
		print OUT "$m\t$start\t$end\t0\n";
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
Usage:    bed2_1k-cov_bdi.pl -i mapped_bed(frombam) -o normalized_rpm_1k -d total_reads -p centromere_position
Function: Calculate mapping reads to 10kb RPM, position unit kb
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d total reads
		  -p position include chr start end chr_length
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-01-26
Notes:    
\n/
    )
}