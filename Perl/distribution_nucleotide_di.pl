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
my %motif;
my $count;
my @di=("AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA","CT","CG","CC");
while (<INPUT>) {
	chomp;
	next if />/;
	$seq=$_;
	$count++;
	foreach  $di(@di) {
		while ($seq=~/$di/g) {
			$pos=pos($seq);
			$motif{$pos}{$di}++;
		}
	}
}
#print OUTPUT "Position\tAA\tAT\tAG\tAC\tTA\tTT\tTG\tTC\tGA\tGT\tGG\tGC\tCA\tCT\tCG\tCC\tAT+TA\tCG+GC\n";
print OUTPUT "Position\tAT+TA\tCG+GC\n";
foreach  $position(sort {$a<=>$b;} keys %motif) {
	print OUTPUT "$position";
#	print $count,"\n";
#	foreach  $di(@di) {
#		printf OUTPUT "\t%.3f",$motif{$position}{$di}/$count;
#	}
	printf OUTPUT "\t%.4f",($motif{$position}{AT}+$motif{$position}{TA})/$count;
	printf OUTPUT "\t%.4f",($motif{$position}{CG}+$motif{$position}{GC})/$count;
	print OUTPUT "\n";

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
Usage:    distribution_nucleotide.pl -i fasta_database -o nucleotide_number_in_position -d nucleotide_list -h header num
Function: match all nucleotide position from batch sequence
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d motif sequence
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-06-03
Notes:    
\n/
    )
}