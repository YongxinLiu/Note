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
$opts{h}=1 unless defined($opts{h});

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
#chrm0    snppos1          ref2     mat_gtyp3        pat_gtyp4        c_gtyp5  phase6   mat_all7 pat_all8 cA9      cC10      cG11      cT12      winning SymCls  SymPval BindingSite     cnv
#1       4648    C       A       C       M       PHASED  C       A       0       11      0       0       M       Asym    0.0009765625    -1      0.902113
open OUTPUT,">$opts{o}";
#chrm    snppos          ref     mat_gtyp        pat_gtyp        c_gtyp  phase   mat_all pat_all cA      cC      cG      cT      winning SymCls  SymPval BindingSite     cnv
#1       4648    C       A       C       M       PHASED  C       A       0       11      0       0       M       Asym    0.0009765625    -1      0.902113

my %count;
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$count{A}=$tmp[9];
	$count{C}=$tmp[10];
	$count{G}=$tmp[11];
	$count{T}=$tmp[12];
	print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$count{$tmp[2]}\t$count{$tmp[3]}\n";
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
Usage:    template.pl -i inpute_file -o output_file -d database	-h header num
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