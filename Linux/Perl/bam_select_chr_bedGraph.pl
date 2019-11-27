#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:r:g:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{r}=1000000 unless defined($opts{r});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my (@tmp1,@tmp2);
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @tmp1,$tmp[0];
}
close DATABASE;



###############################################################################
#Main text.
###############################################################################
#select bam chromosome, and format to bed
my $bed_list;
foreach $chr(@tmp1) {
	`bamtools convert -in $opts{i} -region $chr -format bed -out $chr\.bed`;
	$bed_list.="$chr\.bed ";
}
#merge each chr, and delete temp
`cat $bed_list>$opts{o}\.bed`;
`rm $bed_list`;
$scale=1000000/$opts{r};
`bedtools genomecov -i $opts{o}\.bed -g $opts{d} -scale $scale -bg>$opts{o}\.wig`;
`rm $opts{o}\.bed`;

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
Usage:    bam_select_chr_wig.pl -i inpute_file_bam -o output_file_wig -d database_chr
Function: format bam to wig, and select chr by database list
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -r sample size or mapped reads, use for scale to RPM, default=1000000, no scale
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-05-05
Notes:    Don't run in the same directory at the same time, tmp file will cover
\n/
    )
}