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
$opts{f}=0 unless defined($opts{f});
$opts{s}-=$opts{f};
$opts{e}+=$opts{f};


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#0chr	1strand	2position	3type	4	5level	6mC		7TotalC
#1       C       113447  CG      CG      1.0     7       7
#1       G       113448  CG      CG      1.0     10      10
#1       C       113450  CHG     CT      1.0     7       7
open OUTPUT,">>$opts{o}";
#chr	start	end	C	CG	CHG	CHH
#1	1000	10000	0.90	0.93	0.85	0.02
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my (%C,%mC); # bin total mC level, bin total mC, slid windows mC level, slid windows count
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[0] eq $opts{c}; # not the same chromosome, next
	next if ($tmp[2]<=$opts{s});
	if ($tmp[2]>$opts{e}) {
		printf OUTPUT "$opts{c}\t$opts{s}\t$opts{e}\t%.4f\t%.4f\t%.4f\t%.4f\n",($mC{CG}+$mC{CHG}+$mC{CHH})/($C{CG}+$C{CHG}+$C{CHH}),$mC{CG}/$C{CG},$mC{CHG}/$C{CHG},$mC{CHH}/$C{CHH};
		last;
	}
	$C{$tmp[3]}++;
	$mC{$tmp[3]}+=$tmp[5];
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
Usage:    bin_count_BS1.pl -i CGmap -o mC_average -c chromosome -s start -e end
Function: count assinged region average methyl level, include C, CG, CHG, CHH
Command:  -i inpute file name, *.CGmap (Must)
          -o output file name (Must)
          -h header number, default 0
          -c chromosome
          -s start
          -e end
          -f flank, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-10-22
Notes:    Output C, CG, CHG, CHH level
\n/
    )
}