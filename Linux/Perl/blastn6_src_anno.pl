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
#>ID=SRC1;length_type=24-nt      ath-rep1|Chr1|1|107|-|100|Repbase2012|DNA|ATREP18;
#CCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTAATTCCCTAAACCCGA
my $id;
my %database; #open ath database fasta
while (<DATABASE>) {
	chomp;
	if (/>/) {
		my @tmp=split/\t/;
		$tmp[0]=~s/>//g;
		$id=$tmp[0];
		$database{$tmp[0]}{anno}=$tmp[1];
#		print $chr,"\n";
	}else{
		$database{$id}{len}=length($_);
	}
}
close DATABASE;
my %query; #open osa query fasta
open DATABASE,"<$opts{q}";#reads fasta in hash
while (<DATABASE>) {
	chomp;
	if (/>/) {
		my @tmp=split/\t/;
		$tmp[0]=~s/>//g;
		$id=$tmp[0];
		$query{$tmp[0]}{anno}=$tmp[1];
#		print $chr,"\n";
	}else{
		$query{$id}{len}=length($_);
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
	print OUTPUT "$_\t$query{$tmp[0]}{anno}\t$query{$tmp[0]}{len}\t$database{$tmp[1]}{anno}\t$database{$tmp[1]}{len}\n";
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
Usage:    blastn6_src_anno.pl -i blastn6format -o output.tsv -d genome.fasta
Function: add blast result anno
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file fasta
          -q query file fasta
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-04-15
Notes:    
\n/
    )
}