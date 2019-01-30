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
#Chr1    TAIR10  chromosome      1       30427671        .       .       .       ID=Chr1;Name=Chr1
#Chr1    TAIR10  gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
#Chr1    TAIR10  mRNA    3631    5899    .       +       .       ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
#Chr1    TAIR10  protein 3760    5630    .       +       .       ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
#Chr1    TAIR10  exon    3631    3913    .       +       .       Parent=AT1G01010.1
#Chr1    TAIR10  five_prime_UTR  3631    3759    .       +       .       Parent=AT1G01010.1
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
	if ($tmp[2]=~/gene/i) {
		$tmp[8]=~/Name=([^;]+)/;
		$gene_id=$1;
	}
	if ($tmp[2]=~/exon/i) {
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
Usage:    format_gff_gtf_ath.pl -i inpute_file_gff3 -o output_file_gtf -d database	-h header num
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