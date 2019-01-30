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
#Isolate unique id	Phylum	Familly	Genus	Isolation method	Isolation media	A.thaliana genotype	Soil type	Plant ID
#Root100	Proteobacteria	Phyllobacteriaceae	Aminobacter	Limiting dilution	TSB	Sha	Cologne	plant_202
#Root101	Actinobacteria	Intrasporangiaceae	Phycicoccus	Limiting dilution	TSB	Sha	Cologne	plant_202
my %database; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}=$_;
}

close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#denovo0 Root332
#denovo1 Root217 Root219 Root402 Root568 Root411 Root473 Root70
#denovo2 Root483D1       Root670
	chomp;
	my @tmp=split/\t/;
	my $tmp;
	foreach  (@tmp) {
		if (/denovo/) {
			$tmp=$_;
		}else{
			print OUTPUT "$database{$_}\t$tmp\n";
		}
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
Usage:    pick_otu_anno.pl -i otu_list -o otu_list_taxonomy -d taxonomy	-h header num
Function: format otu to list, add each taxonomy
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-08-12
Notes:    
\n/
    )
}