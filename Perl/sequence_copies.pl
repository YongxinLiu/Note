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
open DATABASE,"<$opts{i}";
my %db; #database in hash
while (<DATABASE>) {
	next if /^>/;
	chomp;
	$db{$_}++;
	$seq=$_;
	$seq=~tr/ACGT/TGCA/;
	$seq=reverse($seq);
	$db{$seq}++;
#	print "$_\t$db{$_}\n";
}

close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#>10_0_20
#GGATTTTTGGAGGATTCGTC
open OUTPUT,">$opts{o}";
#chr0	start1	end2	copies3
#2       20319180        20319200        3
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	$_=~/>([\w_]+)_(\d+)_(\d+)/;
	$seq=<INPUT>;
	chomp($seq);
#	print "$1\t$2\t$3\t$db{$seq}\n";
#	print OUTPUT "$1\t$2\t$3\t$db{$seq}\n";
	$end=$2+1;
	print OUTPUT "$1\t$2\t$end\t$db{$seq}\n";

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
Usage:    sequence_copies.pl -i Position_sequence.fa -o postion_count.bed 
Function: calculate a position fasta into count
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-03-16
Notes:    
\n/
    )
}