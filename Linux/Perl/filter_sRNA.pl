#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:l:r:c:e:n:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{l}="nohup.out" unless defined($opts{l});
$opts{r}=1 unless defined($opts{r});
$opts{c}=20 unless defined($opts{c});
$opts{e}=18 unless defined($opts{e});
$opts{n}=26 unless defined($opts{n});


###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my %database;
#seq0					read1	loci2
#ATCCAGCGCACGGTAGCTT     35      2
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[1];
}

close DATABASE;


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
#seq0					read1	loci2
#ATCCAGCGCACGGTAGCTT     35      2
open OUTPUT,">$opts{o}";
open LOG,">>$opts{l}";
my ($total,$notrRNA);
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$total+=$tmp[1];
	next if defined($database{$tmp[0]});
	next if $tmp[1]<=$opts{r};
	next if $tmp[2]>=$opts{c};
	$len=length($tmp[0]);
	next if ($len<$opts{e} || $len>$opts{n});
	$notrRNA+=$tmp[1];
	print OUTPUT $_,"\n";
}
$per=$notrRNA/$total*100;
#print "Total\tFiltered\tPercentage\n";
$id=basename($opts{i});
printf LOG "$id\t$total\t$notrRNA\t%.2f%%\n",$per;
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
Usage:    filter_sRNA.pl -i sRNA_loci -d trRNA.loci -o sRNA_loci_notrRNA
Function: Filter database contained sRNA from a mapped sRNA_loci file
Command:  -i sample sRNA loci file(Must)
          -o sample sRNA loci file no trRNA (Must)
          -d filtered database
          -r reads count, default > 1
          -c loci cut, default <20
          -e length min cut, default 18
          -n length min cut, default 26
          -l log file
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v2.0
Update:   2013-07-19
Notes:    Add length, read, loci filter
\n/
    )
}