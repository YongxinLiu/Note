#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
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
#>cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
#UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
open INPUT,"<$opts{i}";
my %database;
open OUTPUT,">$opts{o}";
while (<INPUT>) {
        chomp;
        if (/>/) {
            $_=~/>([\w\-\.]+)/;
            $chr=$1;
#           print $chr,"\n";
        }else{
			$_=~s/U/T/g;
            $database{$chr}.=$_;
        }
}
foreach  (keys %database) {
	if (/$opts{s}/) {
		print OUTPUT ">$_\n$database{$_}\n";
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
Usage:    select_miRNA2.pl -i all_miRNA -o species_miRNA
Function: select assign species all miRNA from miRbase
Command:  -i mature.fa (Must)
          -o output file name (Must)
          -s species like ath, osa, gma, zma, mmu
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v2.0
Update:   2013-07-20
Notes:    Output ID only miRNA
\n/
    )
}