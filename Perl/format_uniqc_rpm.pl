#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my @tmp1; #database in array
while (<DATABASE>) {
	chomp;
#	1_9033
#	2_8779
#	3_7624
	my @tmp=split/\t/;
	#print "$tmp[0]\n";
	push @tmp1,$tmp[0];
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#   5820 1_9033
#   1103 2_8779
#    428 3_7624
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my %reads;
while (<INPUT>) {
	chomp;
	my @tmp=split/\s+/;
	#print "$tmp[2]\t$tmp[1]\n";
	$reads{$tmp[2]}=$tmp[1];	
}
foreach (@tmp1) {
#	print $_,"\n";
	$reads{$_}=0 unless defined($reads{$_});
#	printf "$_\t%.2f\n",$reads{$_}/$opts{t}*1000000;
	printf OUTPUT "$_\t%.2f\n",$reads{$_}/$opts{t}*1000000;
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
Usage:    format_uniqc_rpm.pl -i inpute_file -o output_file -d database	-h header num
Function: format unique reads to RPM
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -t total reads for normalization
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-11-03
Notes:    
\n/
    )
}