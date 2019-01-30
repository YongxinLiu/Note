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
print "Database file is $opts{d}\n" if defined($opts{d});

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
open OUTPUT,">$opts{o}";
my %list;
my @filelist=glob "$opts{i}";
foreach $file(@filelist){
	open DATABASE,"<$file";
#	print $file,"\n";
	$file=basename($file);
#	print $file,"\n";
	while (<DATABASE>) {
		my @tmp=split/\t/;
		$list{$file}{nr}++;
		$list{$file}{total}+=$tmp[1];
	}
	close DATABASE;
#	print $list{$file}{total},"\n";
	print OUTPUT $file,"\t",$list{$file}{total},"\t",$list{$file}{nr},"\n";
}
close OUTPUT;

###############################################################################
#Main text.
###############################################################################
#open OUTPUT,">$opts{o}";
#foreach  (sort keys %list) {
#	print OUTPUT "$_\t$list{$_}{total}\t$list{$_}{nr}\n";
#}
#close OUTPUT;

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
Usage:    summary_loci.pl -i inpute_file -o output_file
Function: Summary loci file, and output each loci file contained total reads and non-redundancy reads
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-23
Notes:    
\n/
    )
}