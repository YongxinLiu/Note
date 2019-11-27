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
my %database; #database in hash
while (<DATABASE>) {
#GeneID0		transcriptID1	old geneID, transcript ID is unique
#TUR00001        TUR_T00001      CUFF.1.1
#TUR00002        TUR_T00002      CUFF.2.1
#TUR00003        TUR_T00003      CUFF.3.1
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[1]}=$tmp[0];
}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my %gene;
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#transcriptID0	transcriptID1	evalue2	score3
#TUR_T00003	LOC_Os05g30190.1	1e-20	95.0
#TUR_T00004	LOC_Os07g18510.1	1e-13	73.3
#TUR_T00006	LOC_Os02g22780.1	3e-30	  128
	chomp;
	my @tmp=split/\t/;
	$gi=$database{$tmp[0]};
	$anno=(split(/\./,$tmp[1]))[0];
	if (!defined($gene{$gi}{anno})) {
		$gene{$gi}{anno}=$anno;
		$gene{$gi}{evalue}=$tmp[2];
		$gene{$gi}{score}=$tmp[3];
	}elsif ($tmp[3]>$gene{$gi}{score}){
		$gene{$gi}{anno}=$anno;
		$gene{$gi}{evalue}=$tmp[2];
		$gene{$gi}{score}=$tmp[3];
	}
}
close INPUT;
foreach  (keys %gene) {
	print OUTPUT "$_\t$gene{$_}{anno}\t$gene{$_}{evalue}\t$gene{$_}{score}\n";
}
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
Usage:    format_cDNA2gene.pl -i cDNA_list -o gene_list -d gene_cDNA_list -h header num
Function: format cDNA ID to gene ID, by gene-cDNA list
Command:  -i cDNA related cDNA list, evalue and scores
          -o gene related gene list, evalue and scores
          -d gene cDNA list
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-11-26
Notes:    
\n/
    )
}