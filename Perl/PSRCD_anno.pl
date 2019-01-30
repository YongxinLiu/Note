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
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
my ($left_max,$right_max,$middle_max)=(0,0,0);
open DATABASE,"<$opts{d}";
#database in hash
my %database;
while (<DATABASE>) { #flank file
#Chr1    omicslab        sRNA_cluster    8       1025    57.05   -       .       ID=SRC1;length_type=24-nt       ;       AT1G01010|Chr1|3631|5899|+|2605|TAIR10|protein_coding_gene|NAC domain containing protein 1;ath-rep2|Chr1|1064|1098|+|38|Repbase2012|Simple_repeat|(CACCCCC)n
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[8]}{left}=$tmp[9];
	$left_max=length($tmp[9]) if length($tmp[9])>$left_max;
	$database{$tmp[8]}{right}=$tmp[10];
	$right_max=length($tmp[10]) if length($tmp[10])>$right_max;
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+);length_type=([^;]+)/;
#	print "$1\t$2\n";
	$id=$1;
	$tmp=$2;
#	$id=~s/SRC/$opts{s}\-SRC/;
	$middle_max=length($tmp[9]) if length($tmp[9])>$middle_max;
	#ID/chr/start/end/strand/value/length/left/middle/right
	print OUTPUT "$id\t$tmp[0]\t$tmp[3]\t$tmp[4]\t$tmp[6]\t$tmp[5]\t$tmp\t$database{$tmp[8]}{left}\t$tmp[9]\t$database{$tmp[8]}{right}\n";
}
print "left_max\tmiddle_max\tright_max\n";
print "$left_max\t$middle_max\t$right_max\n";
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
Usage:    PSRCD_anno.pl -i overlap_anno -d flank_anno -o output_file_for_PSRCD_anno
Function: merge SRC gff file overlap and flank anno
Command:  -i inpute overlap name (Must)
          -o output file name (Must)
          -d flank file name
          -s species, ath, osa, gma, zma
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-04-16
Notes:    Delete species function, species formated in result to gff stage
\n/
    )
}