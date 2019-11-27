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
#id0	chr1 strand2	start3	end4 	length5	  loci6	reads7	reads_nrom8	RPM9		RPKM10	 20-nt	21-nt	22-nt	23-nt    24-nt	propotion
#1       Chr1    -       1       1025    1025    1954    73723   73194.29        50.03   48.81   495.00  3524.00 1953.00 13661.00        52468.60        0.63
open OUTPUT,">$opts{o}";
#Chr0	source1	2type	start3	end4		5score(RPM)	strand6	7CDS_phased	detail8(id,name,biotype)
#Chr1    omicslab        SRC     1       1025    50.03   -       .       ID=ath-SRC1;length_type=24-nt
my $length_type;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if ($tmp[11]/$tmp[8]>=0.5) {
		$length_type="20-nt";
	}elsif ($tmp[12]/$tmp[8]>=0.5) {
		$length_type="21-nt";
	}elsif ($tmp[13]/$tmp[8]>=0.5) {
		$length_type="22-nt";
	}elsif ($tmp[14]/$tmp[8]>=0.5) {
		$length_type="23-nt";
	}elsif ($tmp[15]/$tmp[8]>=0.5) {
		$length_type="24-nt";
	}else{
		$length_type="other";
	}	

	print OUTPUT "$tmp[1]\tomicslab\tSRC\t$tmp[3]\t$tmp[4]\t$tmp[9]\t$tmp[2]\t.\tID=$opts{s}\-SRC$tmp[0];length_type=$length_type\n";
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
Usage:    format_dc6_gff.pl -i inpute_file -o output_file -s species
Function: format defined_cluster version6 result to gff3
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -s species, such as ath,osa, gma, zma, default=ath
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-05
Notes:    
\n/
    )
}