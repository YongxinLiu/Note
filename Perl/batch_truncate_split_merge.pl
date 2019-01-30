#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:b:w:', \%opts );
&usage unless ( exists $opts{d});
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
# print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{b}=352 unless defined($opts{b});
$work_dir=`pwd` unless defined($opts{w});
chomp($work_dir);
print "Working directory is : $work_dir\n";
###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my (@tmp1); #database in array
# raw sequencing data ID list
#1018_A
#1018_B
#1018_C
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @tmp1,$tmp[0];
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
print "rm $work_dir\/temp\/2_combined_seqs.fna \n";
`rm $work_dir\/temp\/2_combined_seqs.fna`;
foreach (@tmp1) {
	print "$_\n";
	print "truncate_fasta_qual_files.py -f $work_dir\/data\/$_.fasta -q $work_dir\/data\/$_.qual -b $opts{b} -o $work_dir\/temp\/0_truncated >> log.txt \n";
	`truncate_fasta_qual_files.py -f $work_dir\/data\/$_.fasta -q $work_dir\/data\/$_.qual -b $opts{b} -o $work_dir\/temp\/0_truncated >> log.txt`;
	print "split_libraries.py -f $work_dir\/temp\/0_truncated\/$_\_filtered.fasta -q $work_dir\/temp\/0_truncated\/$_\_filtered.qual -m $work_dir\/data\/$_\_mapping.txt -s 25 -e 0 -b 14 -l $opts{b} -o $work_dir\/temp\/1_demultiplexed\/$_ >>log.txt \n";
	`split_libraries.py -f $work_dir\/temp\/0_truncated\/$_\_filtered.fasta -q $work_dir\/temp\/0_truncated\/$_\_filtered.qual -m $work_dir\/data\/$_\_mapping.txt -s 25 -e 0 -b 14 -l $opts{b} -o $work_dir\/temp\/1_demultiplexed\/$_ >>log.txt`;
	print "cat $work_dir\/temp\/1_demultiplexed\/$_\/seqs.fna >> $work_dir\/temp\/2_combined_seqs.fna \n";
	`cat $work_dir\/temp\/1_demultiplexed\/$_\/seqs.fna >> $work_dir\/temp\/2_combined_seqs.fna`;
}
print "sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;\$/;/g' $work_dir\/temp/2_combined_seqs.fna > $work_dir\/temp/2_combined_seqs_usearch.fasta \n";
`sed 's/ .*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;\$/;/g' $work_dir\/temp\/2_combined_seqs.fna > $work_dir\/temp\/2_combined_seqs_usearch.fasta`;

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
Usage:    batch_truncate_split_merge.pl -d sequencing_ID_list
Function: batch truncate 16S, split into samples and merge into one file
Command:  -d sequencing ID list
          -b truncate length, default 352
          -w work_dir
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-07-28
Notes:    
\n/
    )
}