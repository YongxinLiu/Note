#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
# &usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
$opts{h}=1 unless defined($opts{h});
$opts{i}='SraRunTable.txt' unless defined($opts{i});
$opts{o}='./' unless defined($opts{o});
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

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
#BioSample_s	BioSampleModel_s	Cultivar_s	Experiment_s	ID_s	InsertSize_l	LibraryLayout_s	Library_Name_s	LoadDate_s	MBases_l	MBytes_l	Plant_s	ReleaseDate_s	Rice_Field_s	Run_s	SRA_Sample_s	Sample_Name_s	collection_date_s	env_biome_s	env_feature_s	env_material_s	geo_loc_name_s	host_s	lat_lon_s	Assay_Type_s	AssemblyName_s	BioProject_s	Center_Name_s	Consent_s	LibrarySelection_s	LibrarySource_s	Organism_s	Platform_s	SRA_Study_s	g1k_analysis_group_s	g1k_pop_code_s	source_s
#SAMN02928825	MIGS/MIMS/MIMARKS.host-associated	<not provided>	SRX659832	<not provided>	250	PAIRED	FieldRun1	2016/6/9	30	20	1	2014/7/26	Ditaler_18	SRR1523634	SRS664834	JE.Lundberg.1	2013/7/8	Rice Paddy	Rice Paddy	Rhizosphere	USA: California, Central Valley	Rice	39.41028 N 121.74667 W	AMPLICON	<not provided>	PRJNA255789	UC DAVIS	public	PCR	METAGENOMIC	root metagenome	ILLUMINA	SRP044745	<not provided>	<not provided>	<not provided>
# 2Cultivar_s品种 5InsertSize_l测序长度	6LibraryLayout_s单或双端 11Plant_s植株编号 13Rice_Field_s田间地块 14Run_s	16Sample_Name_s
open INPUT,"<$opts{i}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print "fastq-dump -A $tmp[14] --split-3 -O $opts{o} --gzip\n";
	`fastq-dump -A $tmp[14] --split-3 -O $opts{o} --gzip`;
	print "rename 's/$tmp[14]/$tmp[16]/g' $tmp[14]*\n";
	`rename 's/$tmp[14]/$tmp[16]/g' $tmp[14]*`;
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
        qq/
Usage:    download_sra.pl -i database -o output_file -d database -h header num
Function: download sra file based on SraRunTable.txt, default 14 is SRRID and 16 is sample_ID
Command:  -i inpute file name, default SraRunTable.txt (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-10-24
Notes:    
\n/
    )
}