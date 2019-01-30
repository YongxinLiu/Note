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
$opts{i}="result" unless defined($opts{i});
$opts{o}="result" unless defined($opts{o});
$opts{d}="gene.list" unless defined($opts{d});
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
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
open INPUT,"<$opts{d}";
#ColCvi24hRefvsColOvuleRef_enriched
#ColCvi24hRefvsColOvuleRef_depleted
my %count;
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
#	print "Bak source files\n";
#	`cp $opts{i}/gene_$tmp[0].txt $opts{i}/gene_$tmp[0].bak`;
#	print "cp $opts{i}/gene_$tmp[0].txt $opts{i}/gene_$tmp[0].bak\n\n";

# R output not include row.names, not need this step
#	print "Correction table header\n";
#	`sed -i ':a;N;\$!ba;s/^/ID\\t/' $opts{i}/gene_$tmp[0].txt`;
#	print "sed -i ':a;N;\$!ba;s/^/ID\\t/' $opts{i}/gene_$tmp[0].txt\n\n";

	print "Add signal peptide anno\n";
	`add_anno.sh -i $opts{i}/gene_$tmp[0].txt -d /mnt/bai/public/ref/ensembl/ath/tair10.sp -o temp/temp`;
	print "add_anno.sh -i $opts{i}/gene_$tmp[0].txt -d /mnt/bai/public/ref/ensembl/ath/tair10.sp.raw -o temp/temp\n";
	`sed -i '1 s/\$/SignalP/' temp/temp`;
	print "sed -i '1 s/\$/SignalP/' temp/temp\n";

	print "Add gene description\n";
#	`add_anno.sh -i temp/temp -d /mnt/bai/public/ref/ath/tair10/tair_gene_description13.txt -o $opts{o}/gene_$tmp[0].txt`;
#	print "add_anno.sh -i temp/temp -d /mnt/bai/public/ref/ath/tair10/tair_gene_description13.txt -o $opts{o}/gene_$tmp[0].txt\n\n\n";
#	修改tair_gene_description13为tair_description，原来的只有2.7万个注释，新版有3.2万，但没有之前的详细
	`add_anno.sh -i temp/temp -d /mnt/bai/public/ref/ath/tair10/tair_description.txt -o $opts{o}/gene_$tmp[0].xls`;
	print "add_anno.sh -i temp/temp -d /mnt/bai/public/ref/ath/tair10/tair_description.txt -o $opts{o}/gene_$tmp[0].xls\n\n";
	`sed -i "1 s/\$/Description/" $opts{o}/gene_$tmp[0].xls`;


}
close INPUT;

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
        qq!
Usage:    annotate_gene.pl -i inpute_file -o output_file -d database -h header num
Function: annotate gene by call add_anno.sh from ath gene annotation
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/6/2
Notes:    
\n!
    )
}