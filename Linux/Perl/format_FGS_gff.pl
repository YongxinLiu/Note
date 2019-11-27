#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in hash
#ID0	PFAM1	Panther2	KOG3	KEGG ec4	KEGG Orthology5	TAIR10 hit name6	TAIR10 hit symbol7	TAIR10 hit defline8	rice hit name9	rice hit symbol10	rice hit defline11	
#AC147602.5_FGT004	PF00316	PTHR11556		3.1.3.37	K01100	AT3G55800.1	SBPASE	sedoheptulose-bisphosphatase	LOC_Os04g16680.1		fructose-1,6-bisphosphatase, putative, expressed	
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	my $id=(split/_T/,$tmp[0])[0];
	$id=~s/FGT/FG/;
	next if defined($database{$id});
	if ($#tmp>4 && $tmp[6] ne "") {
		$tmp[8]="NA" if !defined($tmp[8]) || $tmp[8] eq "";
		$tmp[6]=~s/\.\d//;
		$database{$id}="$tmp[6]&$tmp[8]";
		#print "$id\t$database{$id}\n";
	}elsif ($#tmp>7 && $tmp[9] ne ""){
		$tmp[11]="NA" if !defined($tmp[9]) ||$tmp[11] eq "";
		$tmp[9]=~s/\.\d//;
		$database{$id}="$tmp[9]&$tmp[11]";
		#print "$id\t$database{$id}\n";
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#9       ensembl gene    706080  707213  .       +       .       ID=GRMZM2G415022;Name=GRMZM2G415022;biotype=protein_coding
#9       ensembl gene    709107  712584  .       +       .       ID=AC195959.2_FG004;Name=AC195959.2_FG004;biotype=protein_coding
#Mt      ensembl intron  162746  164587  .       -       .       Parent=GRMZM5G861791_T01;Name=intron.929911
#Mt      ensembl intron  204863  205732  .       +       .       Parent=GRMZM5G834128_T01;Name=intron.929934
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (/gene/) {
		$tmp[8]=~/Name=([\w\.\_]+);biotype=([\w\_]+)/;
		for $i(0..7) {
			print OUTPUT "$tmp[$i]\t";
		}
		$database{$1}="NA" unless defined($database{$1});
		print OUTPUT "ID=$1;Note=$2;Annotation=$database{$1}\n";
	}else{
		print $_,"\n";
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
Usage:    format_FGS_gff.pl -i FGS.gff -d phytozome.annotation -o standard.gff
Function: format maize filter gene set (FGS) to standard gff file, only treatment gene and intron
Command:  -i FGS.gff
          -o standard.gff
          -d phytozome homolog annotation
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-22
Notes:    
\n/
    )
}