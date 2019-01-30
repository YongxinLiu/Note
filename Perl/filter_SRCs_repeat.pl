#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:r:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{r}="0.5" unless defined($opts{r});
print "Filter repeat overlap > $opts{r}\n" if defined($opts{r});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

##database in array
#my (@tmp1,@tmp2);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;


#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	$over_len=0;
	chomp;
	my @tmp=split/\t/;
#1-14 basic, 15left, 16overlap, 17right
#ath-SRC1        Chr1    1       201     +       9.06    0.60    1.04    0.00    0.07    0.03    0.25    0.65    24-nt   ;       ath-rep1|Chr1|1|107|-|107|Repbase13|DNA|ATREP18;        AT1G01010|Chr1|3631|5899|+|3429|TAIR10|protein_coding_gene|NAC domain containing protein 1;ath-rep2|Chr1|1064|1098|+|862|Repbase13|Simple_repeat|(CACCCCC)n
	@overlap=split(";",$tmp[15]);
	foreach  $overlap(@overlap) {
		if ($overlap=~/Repbase/) {
#			print $overlap,"\n";
			@tmp1=split('\|',$overlap);
#			print $tmp1[5],"\n";
			$over_len+=$tmp1[5];
		}
	}
	$rate=$over_len/($tmp[3]-$tmp[2]+1);
#	print $_,"\t",$rate,"\n";
	if ($rate<$opts{r}){
		print OUTPUT $_,"\t",$rate,"\n" ;
#		print $rate,"\n";
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
Usage:    filter_SRCs_repeat.pl -i inpute_file -o output_file
Function: filter SRCs enriched in repeat region
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -r rate of overlap, default is 0.5
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-07
Notes:    
\n/
    )
}