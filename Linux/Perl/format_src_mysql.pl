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

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#id0    chr1 strand2    start3  end4    length5 loci6   reads7  reads_nrom8     RPM9    RPKM10	20nt11			21nt12			22nt13		23nt14		24nt15    propotion16(normalized reads in dominant strand)
#520     Chr1    -       6926810 6926897 88      136     237163  237144.83       162.11  1842.11 80967.00        20663.83        24033.00        5758.00 9939.00 0.89
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	
	$id="$opts{s}\-SRC$tmp[0]";
	for (11..15) {
		$ratio[$_]=$tmp[$_]/$tmp[8];
	}
	$copies=$tmp[7]/$tmp[8];
	#ID/chr/start/end/strand/value /strand propotion /copies /20/21/22/23/24propotion 
	printf OUTPUT "$id\t$tmp[1]\t$tmp[3]\t$tmp[4]\t$tmp[2]\t$tmp[9]\t$tmp[16]\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",$copies,$ratio[11],$ratio[12],$ratio[13],$ratio[14],$ratio[15];
	#continue with anno file length/left/middle/right
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
Usage:    format_src_mysql.pl -i inpute_file -o output_file -s ath
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -s species
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-04-02
Notes:    
\n/
    )
}