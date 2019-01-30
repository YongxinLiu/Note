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
#my @filelist=glob "$opts{i}";
#foreach (@filelist){
#	print "$_\n";
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#Chr0	1		2gene	start3			end4		5		strand6		7		type8(id,name,biotype)
#Chr1    TAIR10  CDS     33992   34327   .       -       0       Parent=AT1G01060.1,AT1G01060.1-Protein;
#Chr1    TAIR10  gene    38752   40944   .       -       .       ID=AT1G01070;Note=protein_coding_gene;Annotation=nodulin MtN21 /EamA-like transporter family protein
open OUTPUT,">$opts{o}";
#ID0	Chr1	start2	end3 strand4	source5	type6 description7
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$rpm=$tmp[5]*($tmp[4]-$tmp[3]+1)/1000;
	printf OUTPUT "$_\t%.2f\n",$rpm;
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
Usage:    template.pl -i inpute_file -o output_file
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-04-02
Notes:    
\n/
    )
}