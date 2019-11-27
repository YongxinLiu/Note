#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:t:l:r:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{t}=0 unless defined($opts{t});
$opts{l}=17 unless defined($opts{l});
$opts{r}=27 unless defined($opts{r});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}
#close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#AAAAAAAAAAAAAAA 7151
open OUTPUT,">$opts{o}";
#>1_7151
#AAAAAAAAAAAAAAA
$i=0;
my $len;
while (<INPUT>) {
	$_=~/(\w+)\t(\d+)/;
	$len=length($1);
	if ($len>$opts{l} && $len <$opts{r} && $2>$opts{t})  {
		$i++;
		print OUTPUT ">$i","_$2\n$1\n";#select length 18-30nt, count >2
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
Usage:    format_sRNA_fa.pl -i sRNA_reads_count -o fasta
Function: Format sRNA reads and count to fasta
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -t count filter, default is >0
		  -l length filter, default is >17bp
		  -r length filter, default is <27bp
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-12-05
Notes:    
\n/
    )
}
