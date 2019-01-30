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
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
# miRNA. Mature sequences have type "miRNA".
#
#chr1	.	miRNA_primary_transcript	12425986	12426106	.	+	.	ID=MI0021869_1;accession_number=MI0021869;Name=mmu-mir-6341
#chr1	.	miRNA	12426016	12426038	.	+	.	ID=MIMAT0025084_1;accession_number=MIMAT0025084;Name=mmu-miR-6341;derives_from=MI0021869_1
#based on pre-miRNA ID is same with mature miRNA derives_from
open OUTPUT,">$opts{o}";
#chrom0 Start1 End2 name3  score4 strand5
#chr22  1000   5000 cloneA   960      + 
my %pre;
my ($start,$end);
while (<INPUT>) {
	next if /^\#/;
	chomp;
	my @tmp=split/\t/;
	if ($tmp[2]=~/primary/) {
		$tmp[8]=~/ID=([\w_]+);Alias=([\w]+);Name=([\w-]+)/;
		print "$1\t$3\n";
		$pre{$1}{id}=$3;
		$pre{$1}{start}=$tmp[3];
		$pre{$1}{end}=$tmp[4];
	}elsif ($tmp[8]=~/Derives/i){
		$tmp[8]=~/Name=([\w\.-]+);Derives_from=([\w_]+)/i;
		$preid=$2;
		print "$1\t$2\n";
		if ($tmp[6] eq '+') {
			$start=$tmp[3]-$pre{$preid}{start}+1;
			$end=$tmp[4]-$pre{$preid}{start}+1;
		}else{
			$start=$pre{$preid}{end}-$tmp[4]+1;
			$end=$pre{$preid}{end}-$tmp[3]+1;
		}
		print OUTPUT "$pre{$preid}{id}\t$start\t$end\t$1\t.\t\+\n";
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
Usage:    format_mirChr_mirPre20.pl.pl -i miR gff3 -o mir Bed
Function: Format mature miRNA chromosome position to pre-miRNA position, based on miRbase gff3 file of miRbase20
Command:  -i miRbase_gff3 (Must)
          -o pre-miRNA as chr, BED (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-22
Notes:    
\n/
    )
}