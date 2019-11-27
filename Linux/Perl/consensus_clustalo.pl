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
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my %seq;
while (<INPUT>) {
	chomp;
	next if />/;
	@tmp=split(//,$_);
	for $i(0..$#tmp) {
		$seq{$i}{$tmp[$i]}++;
	}
}

foreach $i(0..$#tmp) {
	$seq{$i}{A}=0 unless defined($seq{$i}{A});
	$seq{$i}{G}=0 unless defined($seq{$i}{G});
	$seq{$i}{C}=0 unless defined($seq{$i}{C});
	$seq{$i}{T}=0 unless defined($seq{$i}{T});
	$seq{$i}{"-"}=0 unless defined($seq{$i}{"-"});
	if ($seq{$i}{A}>=$seq{$i}{G} && $seq{$i}{A}>=$seq{$i}{C} && $seq{$i}{A}>=$seq{$i}{T} && $seq{$i}{A}>=$seq{$i}{"-"}) {
		printf OUTPUT "A\t%.2f\n",$seq{$i}{A}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T}+$seq{$i}{"-"})*100;
		print "A";
		#print "A\t$seq{$i}{A}\t",($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100,"\t",$seq{$i}{A}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100,"\n";
		next;
	}elsif ($seq{$i}{G}>=$seq{$i}{A} && $seq{$i}{G}>=$seq{$i}{C} && $seq{$i}{G}>=$seq{$i}{T} && $seq{$i}{G}>=$seq{$i}{"-"}) {
		printf OUTPUT "G\t%.2f\n",$seq{$i}{G}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T}+$seq{$i}{"-"})*100;
		print "G";
		#printf "G\t%.2f\n",$seq{$i}{G}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		next;
	}elsif ($seq{$i}{C}>=$seq{$i}{A} && $seq{$i}{C}>=$seq{$i}{G} && $seq{$i}{C}>=$seq{$i}{T} && $seq{$i}{C}>=$seq{$i}{"-"}) {
		printf OUTPUT "C\t%.2f\n",$seq{$i}{C}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T}+$seq{$i}{"-"})*100;
		print "C";
		#printf "C\t%.2f\n",$seq{$i}{C}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		next;
	}elsif ($seq{$i}{T}>=$seq{$i}{A} && $seq{$i}{T}>=$seq{$i}{C} && $seq{$i}{T}>=$seq{$i}{G} && $seq{$i}{T}>=$seq{$i}{"-"}) {
		printf OUTPUT "T\t%.2f\n",$seq{$i}{T}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T}+$seq{$i}{"-"})*100;
		print "T";
		#printf "T\t%.2f\n",$seq{$i}{T}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		next;
	}elsif ($seq{$i}{"-"}>=$seq{$i}{A} && $seq{$i}{"-"}>=$seq{$i}{C} && $seq{$i}{"-"}>=$seq{$i}{G} && $seq{$i}{"-"}>=$seq{$i}{T}) {
		printf OUTPUT "-\t%.2f\n",$seq{$i}{"-"}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T}+$seq{$i}{"-"})*100;
		print "-";
		#printf "T\t%.2f\n",$seq{$i}{T}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		next;
	}
}
print "\n";


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
Usage:    consensus_clustalo.pl -i clustalo_rep1A.txt -o clustalo_rep1A.stat
Function: calculate consenseu percentage in clusteO alignment sequence
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-12-01
Notes:    
\n/
    )
}