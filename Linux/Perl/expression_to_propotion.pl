#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Bio::SeqIO;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################
#open DATABASE,"<$opts{d}";
#my @database;
#while (<DATABASE>) {
#	chomp;
#	push @database,$_;
#}
#close DATABASE;


open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my @loc; #2D array
my $title=<INPUT>;
print OUTPUT $title;
my $first=<INPUT>;
my @tmp=split('\t',$first);
push @loc,[@tmp];
my $geneid=(split('\.',$tmp[0]))[0];
my $lastgeneid=$geneid;
my @total;
my $n=0;
for ($n=2;$n<=7 ;$n++) {
	$total[$n]=$tmp[$n];
}
my $count=0;
while (<INPUT>) {
	@tmp=split('\t',$_);
	$geneid=(split('\.',$tmp[0]))[0];
	if ($geneid eq $lastgeneid) {
		push @loc,[@tmp];
		for ($n=2;$n<=7 ;$n++) {
			$total[$n]=$total[$n]+$tmp[$n];
		}
		$count++;
	}else{
		for (my $m=0;$m<=$count ;$m++) {
			print OUTPUT $loc[$m][0],"\t",$loc[$m][1];
			$n=2;
			while ($n<=7) {
				if ($total[$n]==0) {
					print OUTPUT "\t0";
				}else{
					my $propotion=$loc[$m][$n]*100/$total[$n];
					print OUTPUT "\t",$propotion;
				}
				$n++;
			}
			print OUTPUT "\n";
		}
		@loc=();
		push @loc,[@tmp];
		@total=();
		for ($n=2;$n<=7 ;$n++) {
			$total[$n]=$tmp[$n];
		}
		$lastgeneid=$geneid;
		$count=0;
	}
}

for (my $m=0;$m<=$count ;$m++) {
	print OUTPUT $loc[$m][0],"\t",$loc[$m][1];
	$n=2;
	while ($n<=7) {
		if ($total[$n]==0) {
			print OUTPUT "\t0";
		}else{
			my $propotion=$loc[$m][$n]*100/$total[$n];
			print OUTPUT "\t",$propotion;
		}
		$n++;
	}
	print OUTPUT "\n";
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
Usage:    template.pl -i inpute_file -o output_file
Function: Transfer AS expression to proportion.
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-09-20
Notes:    
\n/
    )
}