#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:l:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{t}=10 unless defined($opts{t}); 

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#SRCs list in database
my %database;
while (<DATABASE>) {
#ID0
#ath-SRC8138
#ath-SRC100
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=1;
}
close DATABASE;
my %reads;
if (defined($opts{l})) {
open DATABASE,"<$opts{l}";
#SRCs list in database
while (<DATABASE>) {
#seq0	reads1	loci2
#TTGACAGAAGATAGAGAGCAC   1897174 3
#TGACAGAAGAGAGTGAGCAC    997586  6
	chomp;
	my @tmp=split/\t/;
	$reads{$tmp[0]}=$tmp[1];
}
close DATABASE;
}


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my %sRNA;
while (<INPUT>) {#record all sRNA high expressed 
	chomp;
#ath-SRC1        ATCCCTAAATACCTAATTCCCTA 1161    ATCCCTAAATACCTAATTCCCTAA        491     ATTTAGGGATTCATGGATGTAGGA        465
	my @tmp=split/\t/;
	next unless defined($database{$tmp[0]});
	$sRNA{$tmp[1]}=$tmp[2];
	$sRNA{$tmp[3]}=$tmp[4];
	$sRNA{$tmp[5]}=$tmp[6];
}
my $i=0;
foreach  (keys %sRNA) {
	if (defined($reads{$_})) {
		next unless $reads{$_}>=$opts{t};
		$i++;
		print OUTPUT ">$i\&$reads{$_}\n$_\n";
		next;
	}else{
		next unless $sRNA{$_}>=$opts{t};
		$i++;
		print OUTPUT ">$i\&$sRNA{$_}\n$_\n";
		next;
	}
}


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
Usage:    SRCs_top_reads.pl -i SRCs_reads_database -d  SRCs_list -l reads_count_database -o output_file.fasta
Function: Get top expressed sRNA reads from SRCs ID list
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -l reads_count_database
		  -t reads threshold
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-22
Notes:    
\n/
    )
}