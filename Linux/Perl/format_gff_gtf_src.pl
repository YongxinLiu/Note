#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:f:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{f}=0 unless defined($opts{f});

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
#1       omicslab        SRC     4261    4364    1.82    +       .       ID=zma-SRC1;length_type=24-nt
open OUTPUT,">$opts{o}";
#9       ensembl gene    68562   68582   .       -       0       gene_id "GRMZM2G354611"; transcript_id "GRMZM2G354611_T01";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my ($gene_id,$transcript_id);
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+)/;
	$gene_id=$1;
#	print "$gene_id\n";
	$tmp[8]=~/length_type=([^;]+)/;
	$transcript_id=$gene_id."_".$1;
	$start=$tmp[3]-$opts{f};
	$start=1 if $tmp[3]-$opts{f}<=0;
	print OUTPUT "$tmp[0]\t$tmp[1]\tgene\t",$start,"\t",$tmp[4]+$opts{f},"\t\.\t$tmp[6]\t\.\tgene_id \"$gene_id\"\; transcript_id \"$transcript_id\"\;\n";
	print OUTPUT "$tmp[0]\t$tmp[1]\texon\t",$start,"\t",$tmp[4]+$opts{f},"\t\.\t$tmp[6]\t\.\tgene_id \"$gene_id\"\; transcript_id \"$transcript_id\"\;\n";
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
Usage:    format_gff_gtf_src.pl -i inpute_file_gff3 -o output_file_gtf -d database	-h header num
Function: format gff to gtf of SRCs, used for cufflink
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -f  flank add, default 0 bp
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-11-10
Notes:    
\n/
    )
}