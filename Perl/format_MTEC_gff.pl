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
#1       MTEC    repeat  3619    3931    .       -       .       class=I;subclass=1;order=LTR;superfamily=[RLX] Unknown;family=loukuv;type=Type I Transposons/LTR;name=RLX_loukuv_AC197842_5799
#1       MTEC    repeat  10872   15310   .       -       .       class=II;subclass=1;order=TIR;superfamily=[DTC] CACTA;type=Type II Transposons/TIR;name=ZM_CACTA_42
open OUTPUT,">$opts{o}";
#1       rep13   repeat  147210  147279  560     -       .       ID=zma-rep89;Note=SINE/tRNA;Annotation=SINE2-1_ZM
my $i=1;
my ($id,$note,$anno);
	while (<INPUT>) {
	chomp;
	$id="zma-mtec$i";
	my @tmp=split/\t/;
	if ($tmp[8]=~/;family/) {
		$tmp[8]=~/superfamily=\[(\w+)\]\s+([\w\-_]+);family=([\w\-_]+);[^;]+;name=([\w\_-]+)/;
		$note="$1\&$2\&$3";
		$anno=$4;
	}else{
		$tmp[8]=~/superfamily=\[(\w+)\]\s+([\w\-_]+);[^;]+;name=([\w_\-]+)/;
		$note="$1\&$2\&NA";
		$anno=$3;
	}
	for (0..7) {
		print OUTPUT "$tmp[$_]\t";
	}
	print OUTPUT "ID=$id;Note=$note;Annotation=$anno\n";
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
        qq/
Usage:    format_MTEC_gff.pl -i MTEC.gff -o standard.gff
Function: format MTEC gff to our standard gff
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-22
Notes:    
\n/
    )
}