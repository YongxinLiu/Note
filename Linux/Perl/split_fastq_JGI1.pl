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
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;

#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#@HISEQ03:291:D1427ACXX:6:2308:18427:200568 1:N:0:TTAGGC
#TGCTCCACCATCTCCAAGGGTGGACAGGGTGAACTAGTCAGTCTATGTAGTGCCTAACATCTTTGCCATAAAACTATACGGATCGGCGGGGTCGTCGTTGTTTTTTTTGGTGGGTTTGTTGGTTCTTTGTGTTGTTTTGGGGTGGTTGGG
#+
#4=@DDD++++4+3222+4=FCGH2+14:1:C10*10*0001900019D09BF##################################################################################################
#@HISEQ03:291:D1427ACXX:6:2308:18427:200568 2:N:0:TTAGGC
#CTGGTTGAAGGTGCTAGCAACACTGCGGTATGTCTTGTGGAGTGGGAATGGGGGGGCTCCGGAGAGGGATGCGGTAGTTGGTGGGGTTGGGGGGGGGTCGGGGCGTGGGGGTGGGGGGGGGGGTGGGGTGGTGTTGGGGGTGGTTGGGGA
#+
#+4+4ADD22323C+4,3233+3+3<+)1:)1CH*:CFDF###############################################################################################################
open OUTPUT1,">$opts{o}_1.fastq";
open OUTPUT2,">$opts{o}_2.fastq";
# output_1.fastq
#@HISEQ03:291:D1427ACXX:6:2308:18427:200568 1:N:0:TTAGGC
#TGCTCCACCATCTCCAAGGGTGGACAGGGTGAACTAGTCAGTCTATGTAGTGCCTAACATCTTTGCCATAAAACTATACGGATCGGCGGGGTCGTCGTTGTTTTTTTTGGTGGGTTTGTTGGTTCTTTGTGTTGTTTTGGGGTGGTTGGG
#+
#4=@DDD++++4+3222+4=FCGH2+14:1:C10*10*0001900019D09BF##################################################################################################
# output_2.fastq
#@HISEQ03:291:D1427ACXX:6:2308:18427:200568 2:N:0:TTAGGC
#CTGGTTGAAGGTGCTAGCAACACTGCGGTATGTCTTGTGGAGTGGGAATGGGGGGGCTCCGGAGAGGGATGCGGTAGTTGGTGGGGTTGGGGGGGGGTCGGGGCGTGGGGGTGGGGGGGGGGGTGGGGTGGTGTTGGGGGTGGTTGGGGA
#+
#+4+4ADD22323C+4,3233+3+3<+)1:)1CH*:CFDF###############################################################################################################
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my $i=0;
while (<INPUT>) {
	$i++;
	if ($i<2) {
		chomp;
		my @tmp=split/\s+/;
		print OUTPUT1 "$tmp[0]\n";
	}elsif ($i<5){
		print OUTPUT1 $_;
	}elsif ($i<6){
		chomp;
		my @tmp=split/\s+/;
		print OUTPUT2 "$tmp[0]\n";
	}elsif ($i<8){
		print OUTPUT2 $_;
	}elsif ($i==8){
		print OUTPUT2 $_;
		$i=0;
	}
	
	
}
close INPUT;
close OUTPUT1;
close OUTPUT2;

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
Usage:    split_fastq_JGI.pl -i pair-end.fastq -o output_1.fastq output_2.fastq
Function: split pair-end in one file, pair1 and pair2 in column, split each 4 line one in file
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-12-21
Notes:    format name to same
\n/
    )
}