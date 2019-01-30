#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{t}=1 unless defined($opts{t});
$i=$opts{t};
###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{i}";
#0			1ID				2 transcriptID		3 protein ID		4PFAM 5Panther 	6KOG	7KEGG ec	8KEGG Orthology	9GO	10TAIR10 hit name	11TAIR10 hit symbo	12TAIR10 hit defline	13rice hit name	14rice hit symbol	15rice hit defline	
#31008110	GRMZM2G019411	GRMZM2G019411_T01	GRMZM2G019411_P01	PF01095	PTHR22931,PTHR22931:SF5		3.1.1.11	K01051	GO:0030599,GO:0042545,GO:0005618	AT3G29090.1	ATPME31,PME31	pectin methylesterase 31	LOC_Os10g26680.1		pectinesterase, putative, expressed
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$id=$tmp[$i];
	#print $i,"\t",$id,"\n";
	next if defined($database{$id});
	if ($#tmp>7 && $tmp[10] ne "") {
		$tmp[12]="NA" if !defined($tmp[10]) || $tmp[12] eq "";
		$database{$id}="$tmp[10]&$tmp[12]";
	}elsif ($#tmp>10 && $tmp[13] ne ""){
		$tmp[15]="NA" if !defined($tmp[13]) ||$tmp[15] eq "";
		$database{$id}="$tmp[13]&$tmp[15]";
		#print "$id\t$database{$id}\n";
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
foreach  (keys %database) {
	print OUTPUT "$_\t$database{$_}\n";
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
Usage:    simple_phytzome_anno.pl -i phytozome.annotation -o gene_descript.txt
Function: format maize filter gene set (FGS) to standard gff file, only treatment gene and intron
Command:  -i phytozome homolog annotation
		  -o gene and description
		  -t ID column, phytozome 10 format gene in 1 and transcript in 2, default=1
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-12-08
Notes:    
\n/
    )
}