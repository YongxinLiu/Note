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
#	$database{$tmp[0]}=$tmp[0];
#}
#
#close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
#0ID		1start 2end	3ID									4score	5chr
#scaffold14933	15	116	HWI-ST833:185:7:1101:5286:126295#0	0	+
#scaffold14933	15	115	HWI-ST833:185:7:1302:19015:92455#0	0	+
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
$pointer="";
$count=0;
my ($start,$len);
while (<INPUT>) {
	chomp;
#0				1			2			3			4	5		6		7		8									9
#scaffold14101   .       repeat_region   6101    6263    .       +       .       Name=SSU-rRNA_Ath;class=rRNA;repeat_consensus=N;type=RNA repeats
#scaffold144868  ensembl rRNA_gene       42      156     .       -       .       ID=gene:EPlTURG00000002642;biotype=rRNA;description=5.8S_rRNA [Source:RNAMMER%3BAcc:5.8S_rRNA];external_name=5.8S_rRNA;logic_name=ncrna_eg
	my @tmp=split/\t/;
	if ($tmp[0] eq $pointer) {
		$count++;
		$len=$tmp[3]-$start+1;
		print OUTPUT "$tmp[0]\t$start\t$tmp[3]\tIGS$count\t$len\t$tmp[6]\n";
		$start=$tmp[4];
	}else{
		$start=$tmp[4];
		$pointer=$tmp[0];
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
Usage:    get_gff_igs.pl -i gff -o filtered_gff -d database	-h header num
Function: Read gff file, and output region between anno, such as IGS between gene
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-07-07
Notes:    
\n/
    )
}