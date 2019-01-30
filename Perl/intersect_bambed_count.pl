#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
my %reads;
open DATABASE,"<$opts{d}";
#chrom0		Start1 End2			name3			score4 strand5	6	7	8		9	10		11	Chr12	start13	 end14	ID15	overlap16
#C167069184	66	166	HWI-ST958:142:7:1103:9969:9897#0/1	8	+	66	166	0,0,0	1	100,	0,	C167069184	74	1162	AA_1	92
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$reads{$tmp[15]}++;
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#scaffold126421  6036    6123    IGS8
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($reads{$tmp[3]})) {
		$rpm=$reads{$tmp[3]}*1000000/$opts{t};
		print OUTPUT "$tmp[3]\t$reads{$tmp[3]}\t$rpm\n";
	}else{
		print OUTPUT "$tmp[3]\t0\t0\n";
	}
}
close OUTPUT;
close INPUT;

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
Usage:    intersect_bambed_count.pl -i bed_list -d intersectbed -o output -t total reads for RPM
Function: get the intersect result into ID, count and rpm, only for fastq mapping result
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d intersectBed 
          -t total reads
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0 count reads and RPM in each region
Update:   2014-06-10
Notes:    
\n/
    )
}
