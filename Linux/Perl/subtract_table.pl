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
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
#ID		count	RPM
#IGS6	550	19.7021048904959
#IGS7	0	0
#IGS8	8	0.286576071134486
my @database;
my $i=0;
while ($opts{h}>0) { #filter header
	<DATABASE>;
	$opts{h}--;
}
while (<DATABASE>) {
	chomp;
	@{$database[$i]}=split/\t/;
	$i++;
}

close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#ID		count	RPM
#IGS6	550	19.7021048904959
#IGS7	0	0
#IGS8	8	0.286576071134486
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
$i=0;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT $tmp[0];
	for $j(1..$#tmp) {
		print "$i\t$j\t$tmp[$j]\t$database[$i][$j]\n";
		printf OUTPUT "\t%.6f",$tmp[$j]-$database[$i][$j];
	}
	print OUTPUT "\n";
	$i++;
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
Usage:    subtract_table.pl -i ChIP -o subtract_table -d Input
Function: substract count
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-10-14
Notes:    
\n/
    )
}