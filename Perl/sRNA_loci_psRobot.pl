#!/usr/bin/perl
use warnings;
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#id&reads		chr	strand		start	end		sequences
#784582&2        1       +       1838    1859    TTGCAGAGTAGAAACATAGACT
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split('\t',$_);
	$database{$tmp[5]}++;
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#unique_reads			reads
#TCGGACCAGGCTTCATTCCCC   336322
open OUTPUT,">$opts{o}";
#unique_reads			reads	loci
#TCGGACCAGGCTTCATTCCCC   336322	103
while (<INPUT>){
	chomp;
	my ($seq,$count)=split('\t',$_);
	if (defined $database{$seq}){
		print OUTPUT "$seq\t$count\t$database{$seq}\n";#seq count loc
	}
}
close INPUT;
close INPUT;
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
Usage:    sRNA_loci_psRobot.pl -i sRNA_reads_count -o sRNA_reads_count_loci -d psRobot_mapping
Function: add sRNA NCBI format file loci number, only perfect mapped selected, loci from mapping-v5 result as database
Command:  -i sRNA_reads_count (Must)
          -o sRNA_reads_count_loci
          -d psRobot_mapping
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.2
Update:   2014-06-06
Notes:    originate name "ath_smRNA_loc.pl"
\n/
    )
}
