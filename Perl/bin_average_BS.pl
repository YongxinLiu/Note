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
$opts{t}="C" unless defined($opts{t});
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
#1       G       113452  CHG     CA      0.9     9       10
#1       G       113453  CHH     CC      0.0     0       10
#1       G       113456  CHH     CA      0.0     0       10
#1       G       113457  CHH     CC      0.0     0       11
#1       C       113458  CG      CG      1.0     6       6
#1       G       113459  CG      CG      1.0     11      11
#1       C       113461  CG      CG      1.0     6       6
#1       G       113462  CG      CG      0.818181818182  9       11
#1       G       113464  CHH     CT      0.0     0       11
open OUTPUT,">$opts{o}";
#position
#1_2000	0.76
#1001_3000 0.82
#2001_4000 0.15
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my ($total,$count,$sum_next,$count_next); # bin total mC level, bin total mC, slid windows mC level, slid windows count
print OUTPUT "\t$opts{t}\n";
my %C;
my %mC;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[0] eq $opts{c}; # not the same chromosome, next
	next if ($tmp[2]<=$opts{s} || $tmp[2]>$opts{e}+$opts{b});
	next unless $tmp[3]=~/$opts{t}/;
	if ($tmp[2]>$opts{s}+$opts{b}) {
		print OUTPUT $opts{s}+1,"_",$opts{s}+$opts{b},"\t",$total/$count,"\t",$mC{CG}/$C{CG},"\t",$mC{CHG}/$C{CHG},"\t",$mC{CHH}/$C{CHH},"\n";
		#print "$total\t$count\t$mC\n";
		$opts{s}+=$opts{b};
		$total=0;
		$count=0;
		%C=();
		%mC=();
	}
	$total+=$tmp[5];
	$count++;
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
Usage:    bin_average_BS.pl -i CGmap -o bin_average -c chromosome -s start -e end -b bin_size -l slid_size
Function: count each slid windows average methyl level
Command:  -i inpute file name, *.CGmap (Must)
          -o output file name (Must)
          -t type, CG, CHG, CHH, all, default is all
          -h header number, default 0
          -c chromosome
          -s start
          -e end
          -b bin size, default 2000
          -l slid size, default 2000, now not support
          -r strand split, default none
          -f flank region, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-10-22
Notes:    
\n/
    )
}