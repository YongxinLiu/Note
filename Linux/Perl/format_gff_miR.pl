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
#Chr0	1		2gene	start3			end4		5		strand6		7		type8(id,name,biotype)
#chr1    .       miRNA   28635   28655   .       +       .       ID=MIMAT0004260_1;accession_number=MIMAT0004260;Name=ath-miR838;derives_from=MI0005394_1
open OUTPUT,">$opts{o}";
#chr0	source1	type2		start3	end4 	5		strand6	7		8
#Chr1    miRbase19  miRNA    78932   79032   .       -       .       ID=ath-miR838;Note=miRNA;Annotation=ath-miR838
while (<INPUT>) {
	next if /^#/;
	chomp;
	my @tmp=split/\t/;
	$tmp[0]=~s/c/C/;
	$tmp[8]=~/Name=([^;]+)/;
	print OUTPUT "$tmp[0]\tmiRbase19\t$tmp[2]\t$tmp[3]\t$tmp[4]\t\.\t$tmp[6]\t\.\tID=$1;Note=$tmp[2];Annotation=$1\n";
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
Usage:    format_gff_miR.pl.pl -i inpute_file -o output_file
Function: format miRbase 19 genome gff3 file to our annotation gff format, select pre and mature miRNA
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-04-10
Notes:    
\n/
    )
}