#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:p:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{p}";
#chr	start	end, unit in kb
#Bd1	37435	38515
my %db; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$db{$tmp[0]}{start}=$tmp[1];
	$db{$tmp[0]}{end}=$tmp[2];
}

close DATABASE;


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
   $site=int($hang[1]/1000);
   $num{$Bd}{$site}++;
 } 
close F1;

open OUT,">$opts{o}";
$site{"Bd1"}=75110;
$site{"Bd2"}=59380;
$site{"Bd3"}=59930;
$site{"Bd4"}=48770;
$site{"Bd5"}=28560;


for($m=1;$m<=5;$m++)
 {$Bd="Bd".$m;
 for($n=0;$n<=$site{$Bd};$n++)
  {$start=$n;
   $end=($n+1);
	next if $n<$db{$Bd}{start} || $n>$db{$Bd}{end};
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
Usage:    bed2_1k-cov_bdi.pl -i mapped_bed(frombam) -o normalized_rpm_1k -d total_reads -p centromere_position
Function: Calculate mapping reads to 10kb RPM, position unit kb
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -p position include chr start end
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-01-26
Notes:    
\n/
    )
}