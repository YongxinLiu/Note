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
open OUTPUT,">>$opts{o}";
my $i=1;
my @data;
my $tmp=<INPUT>;
chomp $tmp;
($data[$i][1],$data[$i][2])=split('\t',$tmp);
print OUTPUT $tmp,"\n";
while (<INPUT>) {
	chomp;
	my $pointer=0;
	$i++;
	($data[$i][1],$data[$i][2])=split('\t',$_);
	for my $a(1..($i-1)) {
		if ($data[$a][2] eq $data[$i][2]) {
			$pointer=1;
		}
	}
	if ($pointer==0) {
		print OUTPUT $_,"\n";
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