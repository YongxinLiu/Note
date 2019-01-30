#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:c:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{c}=1 unless defined($opts{c});
my @filelist;
	if (defined $opts{d}) {
	open INPUT,"<$opts{d}";
	while (<INPUT>){
		chomp;
		$filename=$_;
		$filename=$opts{i}.$filename;
		push @filelist,$filename;
	}

}
###############################################################################
#Main text.
###############################################################################
@filelist=glob "$opts{i}" unless defined($opts{d});
my %seq;#all the sequence and reads of total
my %database;#the sequences and reads of every samples
my ($filename,$id);
my @sample;
#print total reads of each samples
foreach (@filelist){
	open INPUT,"<$_";
	print $_,"\n";
	$filename=basename $_;
	$id=$filename;
	push @sample,$id;
	while (<INPUT>){
		chomp;
		my @tmp=split/\t/;
		$database{$tmp[0]}{$id}=$tmp[$opts{c}];
	}
	close INPUT;
}
#print title
open OUTPUT,">$opts{o}";
print OUTPUT "ID";
foreach (@sample){
	print OUTPUT "\t$_";
}
print OUTPUT "\n";
#print maintext
foreach (keys %database){
	print OUTPUT "$_";
	foreach my $sample (@sample){
		if (defined $database{$_}{$sample}){
			printf OUTPUT "\t%.0f",$database{$_}{$sample};
		}else{
			print OUTPUT "\t0";
		}
	}
	print OUTPUT "\n";
}
close OUTPUT;
print scalar @filelist," files have been treated.\n";

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
Usage:    merge_sRNA.pl -i sRNA_file_list -o sRNA table
Function: Merge sRNAs samples file in ID and value format for psRobot_mir
Command:  -i inpute file name (Must)
          -o output file name
          -d database list, instead input if exist
          -c column select, default 1
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  
Update:   2014-11-18
Notes:    
\n/
    )
}
