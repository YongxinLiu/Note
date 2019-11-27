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
#$opts{r}=1 unless defined($opts{r});
$opts{c}=20 unless defined($opts{c});
$opts{e}=18 unless defined($opts{e});
$opts{n}=26 unless defined($opts{n});


###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
my %database;
#seq0					read1
#ATCCAGCGCACGGTAGCTT     35
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$tmp[1];
}
close DATABASE;

open DATABASE,"<$opts{r}";
#database in hash
my %ncRNA;
#seq0
#ATCCAGCGCACGGTAGCTT
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$ncRNA{$tmp[0]}=0;
	$tmp[0]=reverse $tmp[0];
	$tmp[0]=~tr/ACGT/TGCA/;
	$ncRNA{$tmp[0]}=0;
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#seq0
#ATCCAGCGCACGGTAGCTT 
open OUTPUT,">$opts{o}";
#open LOG,">>$opts{l}";
my ($b,$a);
while (<INPUT>) {
	chomp;
	if (defined($ncRNA{$_})){ # filter t/r/sn/snoRNA
		$b++;
		next;
	}
	if (/N/){ # filter N
		$b++;
		next;
	}
	$a=$_;
	if (defined($database{$a})) {
		print OUTPUT "$a\t$database{$a}\n";
	}
	delete $database{$a};
	$a=reverse $a;
	$a=~tr/ACGT/TGCA/;
	if (defined($database{$a})) {
		print OUTPUT "$a\t$database{$a}\n";
	}
	delete $database{$a};
}
close INPUT;
close OUTPUT;
print $b,"\n";

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
Usage:   filter_sRNA_mapped-rRNA.pl -i mapped -d sRNA_reads_count -r mapped_rRNA -o mapped_notrRNA
Function: Select mapped, but non-rRNA reads
Command:  -i mapped sRNA list
          -o mapped notrRNA reads count
          -d reads count
          -r mapped trsnsnoRNA list
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v2.0
Update:   2014-09-19
Notes:    Add length, read, loci filter
\n/
    )
}