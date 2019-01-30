#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:a:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#Model_name      Type				Short_description
#AT1G01010.1     protein_coding  NAC domain containing protein 1
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	my $id=(split/\./,$tmp[0])[0];
	next if defined($database{$id});
	$database{$id}=$tmp[2];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open ANNO,"<$opts{a}";
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#Chr1    omicslab        sRNA_cluster    44570   44956   26.76   -       .       ID=SRC5;length_type=24-nt
#Chr9	source10 type11	start12	end13	score14	strand15 CDS_phased16	detail17(id,name,biotype)
#Chr1    TAIR10  gene    44677   44787   .       +       .       ID=AT1G01073;Note=protein_coding_gene;Name=AT1G01073    111
open OUTPUT,">$opts{o}";
my %anno;
while (<ANNO>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=(\w+)/;
	my $src_id=$1;
	$anno{$src_id}{dis}=0 unless defined($anno{$src_id}{dis});
	$tmp[17]=~/ID=([^;]+);Note=([^;]+)/;
	$geneid=$1;
	$type=$2;
	#print "$tmp[18]\n$anno{$src_id}{dis}\n";
	if ($tmp[18]>$anno{$src_id}{dis}) {
		$anno{$src_id}{dis}=$tmp[18];
		if ($tmp[6] eq $tmp[15]) {#state
			$anno{$src_id}{state}="overlap";
		}else {
			$anno{$src_id}{state}="antisense";
		}
		$anno{$src_id}{type}=$type;
		$anno{$src_id}{id}=$geneid;
		if (defined($database{$geneid})) {#description
			$anno{$src_id}{des}=$database{$geneid};
		}else{
			$anno{$src_id}{des}="NA";
		}
	}
}
close ANNO;

open INPUT,"<$opts{i}";
while (<INPUT>) {
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#1       omicslab        sRNA_cluster    3758    3861    6.57    -       .       ID=SRC1;length_type=24-nt
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=(\w+)/;
	if (!defined($anno{$1}{state})) {
		$anno{$1}{state}="intergenic";
		print OUTPUT "$_\t$anno{$1}{state}\n";
	}else{
		print OUTPUT "$_\t$anno{$1}{state}\t$anno{$1}{type}\t$anno{$1}{id}\t$anno{$1}{dis}\t$anno{$1}{des}\n";
	}
}
close ANNO;
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
Usage:    intersect_gff_maize1.pl -i gff_cluster -a intersect_gff_result -d gene_description -o result
Function: Show the intersect file as a more comprehensive read file,Reads maize gene annotation from phytozome zma annotation, change transcription ID to gene ID; reads each cluster annotation in array, only keep the biggest overlap; result each cluster in one line.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -a intersect_gff_result
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-08-20
Notes:    
\n/
    )
}