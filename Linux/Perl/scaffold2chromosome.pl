#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:u:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{u}="unknown" unless defined($opts{u});


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
#>C158872691
#TGGTGTGTTAGTTGACACGTTTATAACACCACCTTATATCGGAGAACTTGTTATTTCATGTGACATGATCATGTGCATGGATGCATATTAGTGATTTCCTTGTTGTTGGAAATGAATATACTTGTGATGGTGATCACCCTCATGTGCCTTGCATGCATGCCGGAAAGGAATCCGGGAATGCCTCTGCCTTGGGTAGGAGT
#>C158872731
#AAAAAATTAGGGTTACGGGATAAAAGTCAACCGAGTGCAAAGATCCCAGATCTGGATCTACTGTAGCTGACGAACGCACACGAACCGAGGCGGGTTCGTGCTCGCTGGAGTCAGGCCGGGGCTTCGCCCGAGAGGAAGGGGACGTTCGGATCCGGCTGGCCCGTGGCCGGATCTGGCCGGCGACGAGGTCACGGCGCGCG
my %seq;
while (<INPUT>) {
	chomp;
	if (/>/) {
		$_=~/>(\w+)/;
		$name=$1;
	}
	if (/(^[AGCTN]+$)/) {
		$seq{$name}.=$_;
	}
}
open DATABASE,"<$opts{d}";
my %chr;
my $miss=0;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	if (!defined($seq{$tmp[0]})) {
		print "Unknown scaffold ID!!!$tmp[0]\n";
		$miss++;
	}
	if (!defined($seq{$tmp[1]})) {
		$chr{$tmp[1]}.=$seq{$tmp[0]};
	}else{
		$chr{$tmp[1]}.=$seq{$tmp[0]};
	}
	delete $seq{$tmp[0]};
}
print "Missing $miss scaffold\.\n";
foreach (sort keys %seq) {
	$chr{$opts{u}}.=$seq{$_};
	delete $seq{$_};
	last if length($chr{$opts{u}})>2000000000;
}
foreach (sort keys %seq) {
	$chr{"9A"}.=$seq{$_};
	delete $seq{$_};
}

foreach (sort keys %chr) {
	print OUTPUT ">$_\n$chr{$_}\n";
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
Usage:    scaffold2chromosome.pl -i genome.fa -d scaffold_in_chromosome -o pseudo_molecuar -h header_num
Function: Merge all the scaffold to pseudo molecular chromsome, unknown scaffold in 8
		  AA unknown more than 3G, have problem in samtools, trim to 2G at most
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
		  -u unknown chr name, default unknown
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-10-14
Notes:    
\n/
    )
}