#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################
my @filelist=glob "$opts{i}";
my %seq;#all the sequence and reads of total
my %database;#the sequences and reads of every samples
my ($filename,$id);
my @sample;
#print total reads of each samples
open OUTPUT,">$opts{o}.count";
foreach (@filelist){
	open INPUT,"<$_";
	$filename=basename $_;
	$id=$filename;
	print $id,"\n";
	push @sample,$id;
	print OUTPUT "$id\t";
	my $count=0;
	while (<INPUT>){
		chomp;
		if ($_=~/^([AGCT]{18,26})\t(\d+)/){
			$count+=$2;
			next if $2<2;
			$database{$1}{$id}=$2;
			$seq{$1}+=$2;
		}
	}
	close INPUT;
	print OUTPUT "$count\n";
}
close OUTPUT;
#print title
open OUTPUT,">$opts{o}.title";
print OUTPUT "Sequence";
foreach (@sample){
	print OUTPUT "\t$_";
}
print OUTPUT "\n";
close OUTPUT;
#print maintext
open OUTPUT,">$opts{o}";
foreach (keys %seq){
	next if $seq{$_}<100;
	print OUTPUT "$_";
	foreach my $sample (@sample){
		if (defined $database{$_}{$sample}){
			print OUTPUT "\t$database{$_}{$sample}";
		}else{
			print OUTPUT "\t0";
		}
	}
	print OUTPUT "\n";
}
close OUTPUT;
print scalar @filelist," files have been treated.\n";

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
Usage:    merge_table1.pl -i inpute_file -o output_file
Function: Merge all the non-coding RNA NGS in a sequence-sample database, select reads>1, merge reads>=100 and length between 20 to 24-nt sequences
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0 
          v1.1 change all the rpm value Keep decimal point two effective digital
          v1.2 set sequence length must >17bp, debug \$id can't increase, clear the hash on the end
          v1.2 select samples reads>2, and merge reads>100, length between 17 to 26 bp
Update:   2013-04-01
Notes:    
\n/
    )
}
