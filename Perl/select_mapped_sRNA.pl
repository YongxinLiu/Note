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


###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my %database;
#seq0					read1	loci2
#ATCCAGCGCACGGTAGCTT     35      2
while (<DATABASE>) {
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[2];
}

close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#seq0					read1
#ATCCAGCGCACGGTAGCTT     35
open OUTPUT,">$opts{o}";
my ($total,$unique);
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($database{$tmp[0]})) {
		print OUTPUT "$tmp[0]\t$tmp[1]\t$database{$tmp[0]}";
		$total+=$tmp[1];
		$unique++;
	}
}
print "Total reads $total\nUnique reads $unique\n";
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
Usage:    select_mapped_sRNA.pl -i samplesA_sRNA_reads_count -d all_mapped_sRNA_reads_count_loci -o mapped_samplesA_sRNA_reads_count_loci
Function: Select mapped sRNA, based on total mapped sRNA
Command:  -i sample sRNA (Must)
          -o mapped sample sRNA  (Must)
          -d mapped database
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-01-15
Notes:    
\n/
    )
}