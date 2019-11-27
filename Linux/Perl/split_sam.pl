#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:r:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my %database; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[1];
}

close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
open OUTPUTO,">$opts{r}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#0Qname	1Flag	2Rname	3Pos	4MapQ	5Cigar	6Rnext	7Pnext	8Tlen	9Seq	10Qual
#HWI-ST298:518:D0J4BACXX:8:2107:15998:17454      16      1       761     42      24M     *       0       0       TGACATGGATTATGACTCTTACAA        GGJJJJJJJJJHHHHHFFFFFCCC        AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:24 YT:Z:UU
	my @tmp=split/\t/;
	if (defined($database{$tmp[2]})){
		print OUTPUT $_ ;
	}else{
		print OUTPUTO $_ ;;
	}
}
close INPUT;
close OUTPUT;
close OUTPUTO;

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
Usage:    split_sam.pl -i samfile_no_header -o filtered_sam -d chr_list -h header num
Function: split_sam.pl select database chr mapping result
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name, chr_len, bed
		  -r residence result
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}