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
open INPUT,"<$opts{i}";
#细菌rdp方法 greengene13.8注释结果
#OTU_325 k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Glycomycetaceae;g__Glycomyces;s__harbinensis      0.920
#OTU_324 k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae     0.560
#真菌blast方法 unite7.2注释结果
#OTU_4332        k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__unidentified;g__unidentified;s__unidentified       7e-19   SH107231.07FU_JX981477_reps
#OTU_4331        No blast hit    None    None


open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my $i=0;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[0]=~/(\d+)/;
	$i=$1;
	# 有分类级，但没名字的补unnamed，代表有相似序列，但末命名
	$tmp[1]=~s/k__(|unidentified);/k__Unassigned;/;
	$tmp[1]=~s/p__(|unidentified);/p__Unassigned;/;
	$tmp[1]=~s/c__(|unidentified);/c__Unassigned;/;
	$tmp[1]=~s/o__(|unidentified);/o__Unassigned;/;
	$tmp[1]=~s/f__(|unidentified);/f__Unassigned;/;
	$tmp[1]=~s/g__(|unidentified);/g__Unassigned;/;
	$tmp[1]=~s/s__(|unidentified)$/s__Unassigned/;
	# 无分类级补unknown，代表无相似序列新种属
	$tmp[1]="k__Unassigned" unless ($tmp[1]=~/k__/);
	$tmp[1].=";p__Unassigned" unless ($tmp[1]=~/p__/);
	$tmp[1].=";c__Unassigned" unless ($tmp[1]=~/c__/);
	$tmp[1].=";o__Unassigned" unless ($tmp[1]=~/o__/);
	$tmp[1].=";f__Unassigned" unless ($tmp[1]=~/f__/);
	$tmp[1].=";g__Unassigned" unless ($tmp[1]=~/g__/);
	$tmp[1].=";s__Unassigned" unless ($tmp[1]=~/s__/);
	print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\n";
	$i++;
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
        qq!
Usage:    format_taxonomy2full.pl -i result/rep_seqs_tax_assignments.txt -o result/rep_seqs_tax_assignments.txt.full
Function: format taxonomy to full level
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/9/21
Notes:    Only fit assign_taxonomy.py rdp method result
\n!
    )
}