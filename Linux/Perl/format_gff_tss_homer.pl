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
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#1       ensembl gene    4854    9652    .       -       .       ID=gene:GRMZM2G059865;biotype=protein_coding;description=Uncharacterized protein  [Source:UniProtKB/TrEMBL;Acc:C0P8I2];logic_name=genebuilder
#1       ensembl gene    9882    10387   .       -       .       ID=gene:GRMZM5G888250;biotype=protein_coding;logic_name=genebuilder
#1       ensembl gene    109519  111769  .       -       .       ID=gene:GRMZM2G093344;biotype=protein_coding;logic_name=genebuilder
open OUTPUT,">$opts{o}";
#T3G53400.1     3       19796107        19800107        +       .
#T4G30540.1     4       14921327        14925327        +       .
#T1G20560.1     1       7119817	7123817 -       .
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=gene:([^;]+)/;
	if ($tmp[6] eq '+') {
		print OUTPUT "$1\t$tmp[0]\t$tmp[3]\t$tmp[3]\t$tmp[6]\t\.\n";
	}else{
		print OUTPUT "$1\t$tmp[0]\t$tmp[4]\t$tmp[4]\t$tmp[6]\t\.\n";
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
Usage:    format_gff_tss_homer.pl -i inpute_file_gff3 -o tss -h header num
Function: format gff to tss, used for homer create database of promter
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-08
Notes:    
\n/
    )
}