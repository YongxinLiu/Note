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
open DATABASE,"<$opts{d}";
#>1
#GCACATTCTTATCATTTACATTACTGTCTCCTCTCTCCTGTACAGGTACTTGGATGCAGGTTCCTCCAGGAGGCGGCGGGTAGAAAACGAAAATGGGATGAACTTTATTGAATCTGTAATTTAAAGGGmy %database; #database in hash
while (<DATABASE>) {
	chomp;
	$database{$_}=<DATABASE>;
}

close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
#>0_5120 ClusterID, Number
#AGCT Sequence
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}

my $count;
my $long;

while (<INPUT>) {
	chomp;
	#>Cluster 0
	#0       182nt, >145... at -/95.60%
	#1       176nt, >973... at -/96.02%
	if (/^>/) {
		my @tmp=split/\s+/;
		print OUTPUT ">$id\_$count\n$database{$long}" if defined($id);
		$id=$tmp[1];
		$count=0;
	}else{
		$count++;
		if (/\*/) {
			my @tmp=split/\s+/;
			#print $tmp[2],"\n";
			$long=$tmp[2];
			$long=~s/\.\.\.//;
		}
	}
}

print OUTPUT ">$id\t$count\n$database{$long}";

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
Usage:    cd-hit_cluster_represent.pl -i inpute_file.clstr -o output_file.fa -d merged.fa -h header num
Function: Output cluster representative sequences
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-04-30
Notes:    
\n/
    )
}