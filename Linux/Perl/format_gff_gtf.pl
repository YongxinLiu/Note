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
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

##database in array
#my (@tmp1,@tmp2);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;


#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#9       ensembl gene    66347   68582   .       -       .       ID=GRMZM2G354611;Name=GRMZM2G354611;biotype=protein_coding
#9       ensembl exon    68562   68582   .       -       .       Parent=GRMZM2G354611_T01;Name=GRMZM2G354611_E02
#9       ensembl CDS     68562   68582   .       -       0       Parent=GRMZM2G354611_T01;Name=CDS.10
open OUTPUT,">$opts{o}";
#9       ensembl EXON    68562   68582   .       -       0       gene_id "GRMZM2G354611"; transcript_id "GRMZM2G354611_T01";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my ($gene_id,$transcript_id);
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (/gene/i) {
		$tmp[8]=~/Name=([^;]+)/;
		$gene_id=$1;
	}
	if (/exon/i) {
		$tmp[8]=~/Parent=([^;]+)/;
		$transcript_id=$1;
		for (0..7) {
			print OUTPUT "$tmp[$_]\t";
		}
		print OUTPUT "gene_id \"$gene_id\"\; transcript_id \"$transcript_id\"\;\n";
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
Usage:    format_gff_gtf.pl -i inpute_file_gff3 -o output_file_gtf -d database	-h header num
Function: format gff to gtf, used for cufflink
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-11-10
Notes:    
\n/
    )
}