#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:q:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";#reads fasta in hash
#IDq0	IDd1	similiraty2	coverage3	mismatch4	gap5	qstart6	qend7	dstart8	dend9	E-value10	score11
#IDq	IDd			92.59   27			 2       0       296     322     885     911     0.004   39.9
my %database; #open ath database fasta
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($database{$tmp[0]})) {
		next;
	}else{
		$database{$tmp[0]}{sim}=$tmp[2];
		$database{$tmp[0]}{cov}=$tmp[3];
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#IDq0	IDd1	similiraty2	length3	mismatch4	gap5	qstart6	qend7	dstart8	dend9	E-value10	score11
#IDq	IDd			92.59   27			 2       0       296     322     885     911     0.004   39.9
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
#ID	length	
#TR1     250
	if (defined($database{$tmp[0]})) {
		print OUTPUT "$tmp[0]\t$tmp[1]\t$database{$tmp[0]}{sim}\t$database{$tmp[0]}{cov}\n";
	}else{
		print OUTPUT "$tmp[0]\t$tmp[1]\t0\t0\n";
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
Usage:    blastn6_trinity_anno.pl -i ID_list -o ID_list.anno -d blastn_result
Function: add blast max cov and similirity to gene
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file blastn
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-03-08
Notes:    
\n/
    )
}