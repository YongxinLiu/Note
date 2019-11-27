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
getopts( 'i:o:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################

open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	#Chr1    -       AT1G01030.1     13714   13335   13173   11649
	#Chr1    +       AT1G01040.1     23146   24451   24542   24655 
	my @line=split('\t',$_);
	if ($#line<6) {
		next;
	}
	print OUTPUT "\n$line[0]\t$line[1]\t$line[2]";
	my $n=4;
	while ($n<$#line) {
		if ($line[1] eq "+"){
			print OUTPUT "\t",$line[$n]+1,"\t",$line[$n+1]-1;
		}else{
			print OUTPUT "\t",$line[$n]-1,"\t",$line[$n+1]+1;
		}
		$n=$n+2;
	}
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
Function:
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
