#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:s:e:b:', \%opts );
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
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
#print "start\tend\tbin_number\n";
#print "$opts{e}\t$opts{s}\t$opts{b}\n";
my $bin=($opts{e}-$opts{s})/$opts{b};
#print "$bin\n";
my %count;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next if ($tmp[1]<$opts{s} || $tmp[1]>$opts{e});
#	print "$tmp[0]\t$tmp[1]\n";
	for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
#		print "$tmp[0]\t$i\t",$i+$bin,"\n";
		if ($tmp[1]>=$i && $tmp[1]<$i+$bin) {
			$count{$tmp[0]}{$i}++;
			last;
		}
	}
}

foreach  $id(sort keys %count) {
	print OUTPUT "\t$id";
}
print OUTPUT "\n";
for ($j=$opts{s};$j<$opts{e};$j+=$bin) {
#	print "$j\n$opts{s}\t$opts{e}\t$bin\n";
	print OUTPUT "$j";
	foreach  $id(sort keys %count) {
#		print "$id\t$j\t",$j+$bin,"\n";
		$count{$id}{$j}=0 unless defined($count{$id}{$j});
		print OUTPUT "\t$count{$id}{$j}";
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
Usage:    bin_count.pl -i samples -o count -h header num -s start -e end -b bins
Function: count the number in set scale and split bin
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
          -s start number
          -e end number
          -b bin number
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-10-22
Notes:    
\n/
    )
}