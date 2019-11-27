#!/usr/bin/perl
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use strict;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o});
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{s}=2000 unless defined($opts{s});

###############################################################################
#Main text.
###############################################################################

open INPUT,"<$opts{i}";
#0Chr		1gene		2type	3start	4end	5	6strand	7	8detail
#1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;Name=NAC001;biotype=protein_coding;description=NAC domain-containing protein 1 [Source:UniProtKB/Swiss-Prot%3BAcc:Q0WV96];gene_id=AT1G01010;logic_name=araport11
#1       araport11       mRNA    3631    5899    .       +       .       ID=transcript:AT1G01010.1;Parent=gene:AT1G01010;biotype=protein_coding;transcript_id=AT1G01010.1
open OUTPUT,">$opts{o}";
#chr0			start1  end2              ID3       value4    strand5
#scaffold262259    8       31      TRIUR3_26600        gene/upstream/downstream       -
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless defined($tmp[2]);
	next unless $tmp[2]=~/gene/;
	$tmp[8]=~/gene:([^;]+);\w+=([^;]+)/;
	print OUTPUT "$tmp[0]\t$tmp[3]\t$tmp[4]\t$1\t$2\t$tmp[6]\n";
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
Usage:    format_gff2bed_ensemble_gene.pl -i gff3 file (gene) -o bed file 
Function: input standard gff3 file, and output only gene and flank in bed
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
          -s distance, default = 2000
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-07-14
Notes:    
\n/
    )
}