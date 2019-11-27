#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:s:e:b:c:l:', \%opts );
&usage unless ( exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{b}=20 unless defined($opts{b});
$opts{l}=1 unless defined($opts{l});



###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr length
#1	301476924
open OUTPUT,">$opts{o}";
#chr	start	end	ID	.	strand
#1	0	20	1_0_20	.	+
while (<INPUT>) {
	chomp;
	@tmp=split/\t/;
	$start=0;
	while ($start+$opts{b}<=$tmp[1]) {
		$end=$start+$opts{b};
		print OUTPUT "$tmp[0]\t$start\t$end\t$tmp[0]\_$start\_$end\t\.\t\+\n";
		$start+=$opts{l};
	}

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
Usage:    write_bed_bin1.pl -i genome_chr.len (-c chr -s start -e end) -o output_file.bed -b bin_size -l slid_step
Function: write a constitute bed file from start to end, step by bin
Command:  -i genome_chr.len
          -b bin_size, default=20
          -l slid step, default=1
          -o output file name (Must)
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-03-22
Notes:    bed start from 0, e.g. 3:11-20 should be 3:10-20
\n/
    )
}