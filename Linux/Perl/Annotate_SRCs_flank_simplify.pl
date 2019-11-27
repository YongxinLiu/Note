#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
$opts{d}=2000 unless defined($opts{d});
print "Select flank scale $opts{d}\n" if defined($opts{d});

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
#chr0	1				2		start3	end4	RPM5	strand6	7		ID8								left9	right10
#Chr1    omicslab        SRC     1       201     9.06    +       .       ID=ath-SRC1;length_type=24-nt   ;       AT1G01010|Chr1|3631|5899|+|3429|TAIR10|protein_coding_gene|NAC domain containing protein 1;ath-rep2|Chr1|1064|1098|+|862|Repbase13|Simple_repeat|(CACCCCC)n
open OUTPUT,">$opts{o}";
#ID0				leaf1(diatance+ID+description)			right2
#ath-SRC257      317|AT1G05780|Vacuolar ATPase assembly integral membrane protein VMA21-like domain;1|ath-rep747|A-rich  1433|AT1G05785|Got1/Sft2-like vescicle transport protein family;482|ath-rep748|(AACAA)n
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=(\w+\-\w+)/; #match ID
	$id=$1;
	print OUTPUT "$id\t";
	my @tmp1=split/\;/,$tmp[9]; #split Gene and Repeat
	if (length($tmp1[0])>3){ #exclude blank annotate, analysis gene distance
		my @tmp2=split/\|/,$tmp1[0];
		if ($tmp2[5]<$opts{d}) {
			print OUTPUT "$tmp2[5]\|$tmp2[0]\|$tmp2[8];";
		}
	}
	if (length($tmp1[1])>3){ #exclude blank annotate, analysis distance
		my @tmp2=split/\|/,$tmp1[1];
		if ($tmp2[5]<$opts{d}) {
			print OUTPUT "$tmp2[5]\|$tmp2[0]\|$tmp2[8]";
		}
	}
	print OUTPUT "\t";
	my @tmp1=split/\;/,$tmp[10]; #split Gene and Repeat
	if (length($tmp1[0])>3){ #exclude blank annotate, analysis gene distance
		my @tmp2=split/\|/,$tmp1[0];
		if ($tmp2[5]<$opts{d}) {
			print OUTPUT "$tmp2[5]\|$tmp2[0]\|$tmp2[8];";
		}
	}
	if (length($tmp1[1])>3){ #exclude blank annotate, analysis distance
		my @tmp2=split/\|/,$tmp1[1];
		if ($tmp2[5]<$opts{d}) {
			print OUTPUT "$tmp2[5]\|$tmp2[0]\|$tmp2[8]";
		}
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
Usage:    Annotate_SRCs_flank_simplify.pl -i inpute_file -o output_file
Function: simplify flank annotate region, add scale select
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d scale select, default is 2000bp
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-15
Notes:    
\n/
    )
}