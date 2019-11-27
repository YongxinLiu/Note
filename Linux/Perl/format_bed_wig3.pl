#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:r:m:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{m}=-1 unless defined($opts{m});
$opts{r}=1 unless defined($opts{r});
$opts{s}=1000000 unless defined($opts{s});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my %database;
#sequences0			reads1	loci2
#AAAAAAAAAAAAAACTCGGCC   1       1
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[2];
}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr0	start1	end2			sequences3			reads4		strand5
#Chr1    8       31      AGGTTTAGGGTTTAGGGTTTAGGG        4       -
open OUTPUTP,">$opts{o}+.wig";
open OUTPUTN,">$opts{o}-.wig";
#chr	start	end	expression
my %wig;
my %wig1;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next if $tmp[4]<$opts{r};
	$tmp[4]=$tmp[4]/$database{$tmp[3]};
	if ($tmp[5] eq '+') {
		for ($tmp[1]..$tmp[2]) {
			$wig{$tmp[0]}{$_}+=$tmp[4];
		}
	}elsif ($tmp[5] eq '-'){
		for ($tmp[1]..$tmp[2]) {
			$wig1{$tmp[0]}{$_}+=$tmp[4];
		}
	}
}
foreach  my $key1 (sort keys %wig) {
	foreach my $key2 (sort {$a<=>$b;} keys %{$wig{$key1}}) {
		$wig{$key1}{$key2}=$wig{$key1}{$key2}*1000000/$opts{s};
		printf OUTPUTP "$key1\t$key2\t$key2\t%.2f\n",$wig{$key1}{$key2} if defined($wig{$key1}{$key2});
	}
}
foreach  my $key1 (sort keys %wig1) {
	foreach my $key2 (sort {$a<=>$b;} keys %{$wig1{$key1}}) {
		$wig1{$key1}{$key2}=$wig1{$key1}{$key2}*$opts{m}*1000000/$opts{s};
#		print "$wig1{$key1}{$key2}\n";
		printf OUTPUTN "$key1\t$key2\t$key2\t%.2f\n",$wig1{$key1}{$key2} if defined($wig1{$key1}{$key2});
	}
}
close INPUT;
close OUTPUTP;
close OUTPUTN;

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
Usage:    format_bed_to_wig3.pl -i bed(line4_is_reads) -o output_wig_file -d loci file
Function: result in watson and crick two file, normalize by loci, normalize to RPM
Command:  -i use filter loci bed (Must, default use loci < 20 bed)
          -o output file name (Must)
          -d database file name (Must)
          -r int, default is 1, means filter read < 1 
          -m 1\/-1, default is 1, minus is plus, -1 means minus output in minus
          -s scale, total mapped reads, normalize to RPM, default = 1 M
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v3.0
Update:   2014-06-30
Notes:    
\n/
    )
}