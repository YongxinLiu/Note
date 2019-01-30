#!usr/bin/perl
use strict;
my(@hang,$site,%num,$chr,%site,$n,$start,$end,$m);
%num=();

 open F1,"<$ARGV[0]";   ### bed-format
  while(<F1>)
  {chomp;
   @hang=split;
   $chr=$hang[0];
   $site=int($hang[1]/10000);
   $num{$chr}{$site}++;
  } 
close F1;

open OUT,">>$ARGV[1]";
$site{"chr1"}=30147;
$site{"chr2"}=23791;
$site{"chr3"}=23224;
$site{"chr4"}=24206;
$site{"chr5"}=21795;
$site{"chr6"}=16940;
$site{"chr7"}=17682;
$site{"chr8"}=17537;
$site{"chr9"}=15703;
$site{"chr10"}=14963;

for($m=1;$m<=10;$m++)
 {$chr="chr".$m;
 for($n=0;$n<=$site{$chr};$n++)
  {$start=$n*10000;
   $end=($n+1)*10000-1;
   if(exists $num{$chr}{$n}){print OUT "$chr\t$start\t$end\t$num{$chr}{$n}\n";}
   else {print OUT "$chr\t$start\t$end\t0\n";}
  }
 }
close OUT;
