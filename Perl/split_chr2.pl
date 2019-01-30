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
#Main text.
###############################################################################
open DATABASE,"<$opts{d}";#reads fasta in hash
#>1
#GAATTCCAAAGCCAAAGATTGCATCAGTTCTGCTGCTATTTCCTCCTATCATTCTTTCTG
#ATGTTGAAAATGATATTAAGCCTAGGATTCGTGAATGGGAGAAGGTATTTTTGTTCATGG
my %database;
my $chr;
while (<DATABASE>) {
	chomp;
	if (/>/) {
		$_=~/>(\w+)/;
		$chr=$1;
#		print $chr,"\n";
	}else{
		$database{$chr}.=$_;
	}
}
close DATABASE;

open INPUT,"<$opts{i}";
`mkdir $opts{o}`;
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	if (defined($database{$_})) {
		open OUTPUT,">$opts{o}$_\.fa";
		print OUTPUT ">$_\n$database{$_}\n";
		close OUTPUT;
	}
}

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
Usage:    split_chr2.pl -i chr list -o output_dir\/ -d genome fasta
Function: split fasta genome into each file, only report chrID in database
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d chrID list file name (Must)
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2014-11-20
Notes:    
\n/
    )
}