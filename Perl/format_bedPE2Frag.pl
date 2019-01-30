#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{s}=0 unless defined($opts{s});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr0	start1	end2			ID3				4		strand5
#3       3       103     SRR639498.15149949      1       +
#3       44      144     SRR639498.4065068       1       +
#3       62      162     SRR639498.15149949      1       -
#3       106     206     SRR639498.4065068       1       -
open OUTPUT,">$opts{o}";
#chr0	start1	end2	ID3 width4	strand5
#1	332	2844		SRR639498.15149949	284	+
my %pair;
#print OUTPUT "space\tstart\tend\twidth\tstrand\n";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[3]=(split/\//,$tmp[3])[0];
	if (defined($pair{$tmp[3]}{chr})) {
		next unless $pair{$tmp[3]}{chr} eq $tmp[0];
		print OUTPUT "$tmp[0]\t",$pair{$tmp[3]}{start}-$opts{s},"\t",$tmp[2]-$opts{s},"\t$tmp[3]\t",$tmp[2]-$pair{$tmp[3]}{start}+1,"\t\+\n";
	}else{
		$pair{$tmp[3]}{chr}=$tmp[0];
		$pair{$tmp[3]}{start}=$tmp[1];
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
Usage:    format_bedPE2Frag.pl -i bed(line4_is_reads) -o bed(pair-end as fragment)
Function: format bed(from bam) to nucleR use pair-end as fragment
Command:  -i sorted bed
          -o output file name (Must)
          -d database file name (Opt)
          -s step size, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2015-06-11
Notes:    
\n/
    )
}