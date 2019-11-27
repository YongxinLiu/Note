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

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#Chr	Length
#scaffold9858    1066088
#scaffold7153    874486
#scaffold8955    808725
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
`cp $opts{h} $opts{o}\.sam`;
print "cp $opts{h} $opts{o}\.sam\n";
foreach  (@tmp1) {
	`samtools view $opts{i} "$_">>$opts{o}\.sam\n`;
	print "samtools view $opts{i} \"$_\">>$opts{o}\.sam\n";
}
`samtools view -bS $opts{o}.sam -o $opts{o}.bam`;
print "samtools view -bS $opts{o}.sam -o $opts{o}.bam\n";
`samtools index $opts{o}\.bam`;
print "samtools index $opts{o}\.bam";
`rm $opts{o}.sam`;
print "rm $opts{o}.sam\n";



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
Usage:    bam_select_chr.pl -i inpute_file_bam -o output_file_bam -d scaffold_list -h header
Function: Select some scaffold in bam, use samtools view quickly show each scaffold, need give a new header
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h chromosome number
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-11-07
Notes:    
\n/
    )
}
