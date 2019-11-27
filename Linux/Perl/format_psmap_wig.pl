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

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#id&reads0	chr1	strand2		start3	end4 	sequences5
#21719662&1      1       +       656     677     TGGATATAAACTATTTTTGGCT
open OUTPUTP,">$opts{o}+.wig";
open OUTPUTN,">$opts{o}-.wig";
#chr	start	end	expression

while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	my @tmp1=split/&/,$tmp[0];
	if ($tmp[2] eq '+') {
		for ($tmp[3]..$tmp[4]) {
			$wig{$tmp[1]}{$_}+=$tmp1[1];
		}
	}elsif ($tmp[2] eq '-'){
		for ($tmp[3]..$tmp[4]) {
			$wig1{$tmp[1]}{$_}+=$tmp1[1];
		}
	}
}
foreach  my $key1 (sort keys %wig) {
	foreach my $key2 (sort {$a<=>$b;} keys %{$wig{$key1}}) {
		print OUTPUTP "$key1\t$key2\t$key2\t$wig{$key1}{$key2}\n" if defined($wig{$key1}{$key2});
	}
}
foreach  my $key1 (sort keys %wig1) {
	foreach my $key2 (sort {$a<=>$b;} keys %{$wig1{$key1}}) {
		print OUTPUTN "$key1\t$key2\t$key2\t$wig1{$key1}{$key2}\n" if defined($wig1{$key1}{$key2});
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
Usage:    format_psmap_wig.pl -i psmap_mapping -o wig+\/-
Function: format psRobot mapping result to wig, not normalize by loci
Command:  -i psmap_mapping(Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-03-04
Notes:    
\n/
    )
}