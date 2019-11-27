#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
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
#   3648 3542_1  ata
#   3474 3321_1  ata
#   3442 3849_1  ata
#   3191 3977_1  ata
#   3083 5085_1  ata
#   3045 4682_1  ata
#   2958 3307_1  ata
#   2761 3515_1  ata
#   2667 3579_1  ata
#   2598 4374_1  ata
#   2520 4021_1  ata
#   2441 3482_1  ata
#   2337 4146_1  ata
#   2168 4222_1  ata
#   2137 4980_1  ata
#   2098 3989_1  ata
#   1985 3722_1  ata
#   1885 831_2   tur
#   1860 4512_1  ata
#   1810 1878_1  tur
#   1750 455_1   tur
#   1686 813_1   tur
#   1681 1109_1  tur
#   1605 1064_1  tur
#   1594 1671_1  tur
#   1594 2310_1  tur
#   1584 305_1   tur
#   1578 1299_1  tur

open OUTPUT,">$opts{o}";
# ID	AAcopie	DDcopies	AA/DD
# 1299_1 1578 1578 1
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my %input;
while (<INPUT>) {
	chomp;
	my @tmp=split/\s+/;

	$input{$tmp[2]}{$tmp[3]}=$tmp[1];

}
close INPUT;

foreach  (keys %input) {
	$input{$_}{tur}=0 unless defined($input{$_}{tur});
	$input{$_}{ata}=0 unless defined($input{$_}{ata});
	print OUTPUT "$_\t$input{$_}{tur}\t$input{$_}{ata}\t",($input{$_}{tur}+1)/($input{$_}{ata}+1),"\n";
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
Usage:    format_uniqc_AD.pl -i inpute_file -o output_file -d database	-h header num
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}