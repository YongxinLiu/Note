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

# open a list file
my %list;
my @filelist=glob "$opts{i}";
foreach $file(@filelist){
	# print $file,"\n";
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
	$filename="$_\/fastqjoin.join.fastq"; # 合并目录+/+固定文件名
	if (-e $filename) {
		$_=~/([\w\.]+)_1/; # 目录包括_1
		$id=$1;
		open INPUT,"<$filename";
		my $count=0;
		while (<INPUT>) { #filter header
			$count++;
			$seq=<INPUT>;
			print OUTPUT ">$id\_$count\n$seq";
			<INPUT>;
			<INPUT>;
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
Usage:    merge_fastq-join.pl -i input_dir -o output_file -d samples_list -h header num
Function: merge_fastq-join into one file
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