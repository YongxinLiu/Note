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
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#1050608 L2S309_105322   L2S357_93299    L4S112_137799   L4S112_115484   L4S112_148483   L2S155_106013   L2S382_106113
#2595164 L4S63_140520    L4S63_124712
open OUTPUT,">$opts{o}";
#1050608 L2S357_93299
#1050608 L2S309_105322
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my $count=0;
my %otu;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	for (1..$#tmp) {
#		print @tmp,"\t",$#tmp,"\n";
		print OUTPUT "$tmp[0]\t$tmp[$_]\n";
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
        qq!
Usage:    format_uctxt2list.pl -i inpute_file -o output_file -d database	-h header num
Function: Format OTU wide table to long table. ~/culture/medicago/fungi]$format_uctxt2list.pl -i uclust_100/171225_clean_otus.txt -o uclust_100/171225_clean_otus.list
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n!
    )
}