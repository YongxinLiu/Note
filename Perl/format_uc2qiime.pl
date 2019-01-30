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
#H       23      308     99.7    +       0       0       308M    HN.32.4_0;barcodelabel=HN.32.4; OTU_24
#H       46      260     100.0   +       0       0       260M    HNS.2_1;barcodelabel=HNS.2;     OTU_47
#H       3       273     100.0   +       0       0       273M    LN.F2.8.6_11;barcodelabel=LN.F2.8.6;    OTU_4
open OUTPUT,">$opts{o}";
#337212  L4S112_121780
#1050608 L2S309_105322   L2S357_93299    L4S112_137799   L4S112_115484   L4S112_148483   L2S155_106013   L2S382_106113
#2595164 L4S63_140520    L4S63_124712
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my $count=0;
my %otu;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if ($tmp[9]=~/OTU/){
		@tmp1=split(";",$tmp[8]);
		$count++;
		if (defined($otu{$tmp[9]})) {
			$otu{$tmp[9]}.="\t$tmp1[0]";
		}else{
			$otu{$tmp[9]}=$tmp1[0];
		}

	}
}
foreach  (keys %otu) {
	print OUTPUT "$_\t$otu{$_}\n";
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
Usage:    format_uc_qiime.pl -i inpute_file -o output_file -d database	-h header num
Function: format usearch otu_talbe uc file to qiime used pick_otu txt result
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}