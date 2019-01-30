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
open DATABASE,"<$opts{d}";
#database in hash
my %database;
while (<DATABASE>) {
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#9       ensembl gene    311374  318659  .       -       .       ID=GRMZM2G310569;Note=protein_coding;Annotation=AT3G02850&STELAR K+ outward rectifier
#9       ensembl gene    563756  566860  .       +       .       ID=GRMZM2G341658;Note=protein_coding;Annotation=AT1G14780&MAC/Perforin domain-containing protein
	chomp;
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+);Note=[^=]+=(.*)/;
#	print "$1\t$2\n";
	$database{$1}=$2;
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#>1_105  Score: 2.5      GRMZM2G050485_T01 seq=cdna; coord=8:142225800..142228222:1; parent_gene=GRMZM2G050485
	if (/>/) {
		chomp;
		$_=~/=([\w\.\_]+)$/;
		print $1,"\n$database{$1}\n";
		print OUTPUT $_,"\t$database{$1}\n";
	}else{
		print OUTPUT $_;
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
Usage:    maize_gTP_anno.pl -i input.gTP -d maize phytozome annotation -o annotated.gTP	-h header num
Function: Add maize target file gene description
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-09-10
Notes:    
\n/
    )
}