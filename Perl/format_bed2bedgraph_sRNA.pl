#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:r:m:s:a:b:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{m}=-1 unless defined($opts{m});
$opts{r}=1 unless defined($opts{r});
$opts{s}=1000000 unless defined($opts{s});
$opts{a}=20 unless defined($opts{a});
$opts{b}=25 unless defined($opts{b});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr0	start1	end2	sequences3			reads4		strand5
#IGS-A   10      34      5998225_2       255     -
open OUTPUTP,">$opts{o}";
#chr	start	end	expression
my %wig;
my %wig1;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[3]=~/(\d+)_(\d+)/;
	$tmp[4]=$2;
	$len=$tmp[2]-$tmp[1];
	next if $len<$opts{a};
	next if $len>$opts{b};
	next if $tmp[4]<$opts{r};
	for ($tmp[1]..$tmp[2]) {
			$wig{$tmp[0]}{$_}+=$tmp[4];
	}
}
foreach  my $key1 (sort keys %wig) {
	foreach my $key2 (sort {$a<=>$b;} keys %{$wig{$key1}}) {
		$wig{$key1}{$key2}=$wig{$key1}{$key2}*1000000/$opts{s};
		printf OUTPUTP "$key1\t$key2\t$key2\t%.2f\n",$wig{$key1}{$key2} if defined($wig{$key1}{$key2});
	}
}
close INPUT;
close OUTPUTP;

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
Usage:    format_bed2bedgraph_sRNA.pl -i bam_bed_sRNA -o output_file -s total_reads
Function: result in one file, not normalize by loci, normalized to RPM
Command:  -i use filter loci bed
          -o output file name (Must)
          -d database file name (Opt)
          -a min length, 20 (Opt)
          -b min length, 25 (Opt)
          -s scale, total mapped reads, normalize to RPM, default = 1 M
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2014-03-03
Notes:    
\n/
    )
}