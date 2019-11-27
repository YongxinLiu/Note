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
my %list;
my @filelist=glob "$opts{i}";
foreach $file(@filelist){
	$list{$file}=0;
}

###############################################################################
#Main text.
###############################################################################
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
open OUTPUT,">$opts{o}";
open OUTPUTL,">$opts{o}\.log";
my $id;
foreach (keys %list) {
	chomp;
	$filename=$_; # 合并目录+/+固定文件名
	if (-e $filename) {
		$_=~/([\w\.]+)\.fna/; # 目录包括_1
		$id=$1;
		open INPUT,"<$filename";
		my $count=0;
		while (<INPUT>) { #filter header
			$count++;
			$seq=<INPUT>;
			print OUTPUT ">$id\_$count\n$seq";
		}
		print $id,"\t",$count,"\n";
		print OUTPUTL $id,"\t",$count,"\n";
		
	}else{
		print "Files doesn't exist. $_\/fastqjoin.join.fastq\n";
	}
	
}

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
Usage:    merge_fastq_fasta.pl -i input_dir -o output_file -d samples_list -h header num
Function: merge_fastq into one fasta file
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-11-11
Notes:    
\n/
    )
}
