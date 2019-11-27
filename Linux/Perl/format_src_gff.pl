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
$opts{s}="ath" unless defined($opts{s});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#id0	chr1 		start2	end3 	strand4	RPM5	strand_ratio6	copies7 20-24 length13	overlap14
#ath-SRC2        Chr1    207     354     -       0.87    0.77    1.00    0.01    0.00    0.00    0.07    0.90    24-nt   intergenic
open OUTPUT,">$opts{o}";
#Chr0	source1	2type	start3	end4		5score(RPM)	strand6	7CDS_phased	detail8(id,name,biotype)
#Chr1    omicslab        SRC     1       1025    50.03   -       .       ID=ath-SRC1;length_type=24-nt
my $length_type;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT "$tmp[1]\tomicslab\tSRC\t$tmp[2]\t$tmp[3]\t$tmp[5]\t$tmp[4]\t.\tID=$tmp[0];length_type=$tmp[13]\n";
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
Usage:    format_src_gff.pl -i inpute_file -o output_file -s species
Function: format defined_cluster version6 result to gff3
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -s species, such as ath,osa, gma, zma, default=ath
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-06-25
Notes:    
\n/
    )
}
