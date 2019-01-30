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
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";

open DATABASE,"<$opts{d}";
my @database;
while (<DATABASE>) {
	chomp;
	push @database,$_;
#	print $_,"\n";
}
close DATABASE;

while (<INPUT>) {
	chomp;
	if (/(osa-miR|ath-miR)/) {
		print OUTPUT $_,"\n";
		my @tmp=split('\t',$_);
		my $geneid=(split('\.',$tmp[4]))[0];
			my $n=scalar @database;
			for (my $i=0; $i<$n; $i++) {
				if ($database[$i]=~m/$tmp[0]/) {
					if ($database[$i]=~m/$geneid/) {
						print OUTPUT $database[$i],"\n",$database[$i+2],"\n",$database[$i+3],"\n",$database[$i+4],"\n";
					}
				}
			} 
	}else{
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