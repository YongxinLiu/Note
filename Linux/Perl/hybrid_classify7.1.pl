#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'p:m:h:o:s:', \%opts );
&usage unless ( exists $opts{p} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Paternal file is $opts{p}\nMaternal file is $opts{m}\nHybrid file is $opts{h}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{s}=0.2 unless defined($opts{s});

###############################################################################
#Read the database in memory(opt)
###############################################################################
my (@id,@paternal,@maternal,@hybird);
#id0	loci1	reads2	norm_reads3	rpkm4
#1       8       48      48.00   10.15
open PATERNAL,"<$opts{p}";
while (<PATERNAL>) {
	chomp;
	my @tmp=split/\t/;
	push @id,$tmp[0];
	push @paternal,$tmp[1];
}
close PATERNAL;
open MATERNAL,"<$opts{m}";
while (<MATERNAL>) {
	chomp;
	my @tmp=split/\t/;
	push @maternal,$tmp[1];
}
close HYBRID;
open HYBRID,"<$opts{h}";
while (<HYBRID>) {
	chomp;
	my @tmp=split/\t/;
	push @hybrid,$tmp[1];
}
close HYBRID;

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
#id	paternal	maternal	hybrid	MPV	class
#16090	68.87	0.13	0.53	34.5	MP/HP/LP/AHP/BLP/BMH/BML
my ($mpv,$hp,$lp);
for (0..$#paternal) {
	$mpv=($paternal[$_]+$maternal[$_])/2;
	$hp=($paternal[$_]>$maternal[$_]? $paternal[$_]:$maternal[$_]);
	$lp=($paternal[$_]<$maternal[$_]? $paternal[$_]:$maternal[$_]);
#	print "max=$hp\tmin=$lp\n"
	if ($hybrid[$_]>=($mpv*(1-$opts{s})) && $hybrid[$_]<=($mpv*(1+$opts{s}))) {#MP
		print OUTPUT "$id[$_]\t$paternal[$_]\t$maternal[$_]\t$hybrid[$_]\t$mpv\tMP\n";
	}elsif ($hybrid[$_]>=($lp*(1-$opts{s})) && $hybrid[$_]<=($lp*(1+$opts{s}))) {#LP
		print OUTPUT "$id[$_]\t$paternal[$_]\t$maternal[$_]\t$hybrid[$_]\t$mpv\tLP\n";
	}elsif ($hybrid[$_]>=($hp*(1-$opts{s})) && $hybrid[$_]<=($hp*(1+$opts{s}))) {#HP
		print OUTPUT "$id[$_]\t$paternal[$_]\t$maternal[$_]\t$hybrid[$_]\t$mpv\tHP\n";
	}elsif ($hybrid[$_]<($lp*(1-$opts{s}))) {#BLP
		print OUTPUT "$id[$_]\t$paternal[$_]\t$maternal[$_]\t$hybrid[$_]\t$mpv\tBLP\n";
	}elsif ($hybrid[$_]>($hp*(1+$opts{s}))) {#AHP
		print OUTPUT "$id[$_]\t$paternal[$_]\t$maternal[$_]\t$hybrid[$_]\t$mpv\tAHP\n";
	}elsif ($hybrid[$_]>($lp*(1+$opts{s})) && $hybrid[$_]<($mpv*(1-$opts{s}))) {#BML
		print OUTPUT "$id[$_]\t$paternal[$_]\t$maternal[$_]\t$hybrid[$_]\t$mpv\tBML\n";
	}elsif ($hybrid[$_]>($mpv*(1+$opts{s})) && $hybrid[$_]<($hp*(1-$opts{s}))) {#BMH
		print OUTPUT "$id[$_]\t$paternal[$_]\t$maternal[$_]\t$hybrid[$_]\t$mpv\tBMH\n";
	}
}
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
Usage:    hybrid_classify7.1.pl -p B1_cluster.rpkm -m A2_cluster.rpkm -h A3_cluster.rpkm -o B1_A2_A3.class -s 0.05
Function: Basic on paternal, maternal and hybrid RPM, classify genes, SRCs into 7 groups, middle parent-MP, high parent-HP, low parent-LP, above high parent-AHP, below low parent-BLP, between middle and high-BMH, between middle and low-BML.
Command:  -p paternal RPM
          -m maternal RPM
          -h hybird RPM
          -o output result
		  -s scale for +-$opts{s}, default 1+-0.2
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-08-08
Notes:    Originate from classify_DEC.pl
		  exchange AHP and BLP, HP and LP, AHP and BLP order, check result is change?
\n/
    )
}