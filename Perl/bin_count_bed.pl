#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:c:s:e:b:l:f:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{b}=2000 unless defined($opts{b});
$opts{l}=1000 unless defined($opts{l});
$opts{t}=10000000 unless defined($opts{t});
$opts{f}=0 unless defined($opts{f});
$opts{s}-=$opts{f};
$opts{e}+=$opts{f};


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#0chr	1start	2end	3ID	4score	5strand
#chr1	593	693	HWI-ST298:373:D09UHACXX:6:1102:5407:124774/1	60	+
#chr1	733	833	HWI-ST298:373:D09UHACXX:6:1102:5407:124774/2	60	-
#chr1	809	909	HWI-ST298:373:D09UHACXX:6:1307:15270:180191/1	60	+
open OUTPUT,">$opts{o}";
#position
#1_2000	0.76
#1001_3000 0.82
#2001_4000 0.15
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my ($rpm,$count,$sum_next,$count_next); # bin total mC level, bin total mC, slid windows mC level, slid windows count
print OUTPUT "\tCenH3\n";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[0] eq $opts{c}; # not the same chromosome, next
	$tmp[2]=($tmp[1]+$tmp[2])/2;
	next if ($tmp[2]<=$opts{s} || $tmp[2]>$opts{e}+$opts{b});
	if ($tmp[2]>$opts{s}+$opts{b}) {
#		print "$count\n";
		$rpm=$count/$opts{t}*1000000;
		print OUTPUT $opts{s}+1,"_",$opts{s}+$opts{b},"\t$rpm\n";
		$opts{s}+=$opts{b};
		$count=1;
	}
	$count++;
}


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
Usage:    bin_count_bed.pl -i CGmap -o bin_average -c chromosome -s start -e end -b bin_size -l slid_size
Function: count each slid windows average methyl level
Command:  -i inpute file name, *.CGmap (Must)
          -o output file name (Must)
          -t total reads, default=10000000
          -h header number, default 0
		  -c chromosome
          -s start
          -e end
          -b bin size, default 2000
		  -l slid size, default 1000
		  -r strand split, default none
		  -f flank region, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-04
Notes:    
\n/
    )
}