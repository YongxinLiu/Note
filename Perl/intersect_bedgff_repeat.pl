#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open LOCI,"<$opts{t}";
#sequences0		reads1	loci2
#AAAAAAAAAAAAAACTCGGCC   1       1
my (%loci,$total);
while (<LOCI>) {
	chomp;
	my @tmp=split/\t/;
	$loci{$tmp[0]}=$tmp[2];#preserve all the loci number in hash
	$total+=$tmp[1];
}
close LOCI;
#open DATABASE,"<$opts{d}";
##chrom0	source1	type2 start3  end4		5 strand6		7		ID8															
##9       ensembl gene    66347   68582   .       -       .       ID=GRMZM2G354611;Name=GRMZM2G354611;biotype=protein_coding
#my (%id,%len);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$tmp[8]=~/Name=([^;]+)/;
##	print $1,"\n";
#	$id{$1}=$1;
#	$len{$1}=$tmp[4]-$tmp[3]+1;
#}
#close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
my (%reads_norm,,%rpm,%rpkm,%length);
my ($start,$end)=(0,0);
my $family;
while (<INPUT>) {
#chrom0 start1 end2				name3			reads4 strand5 chrom6	source7	type8 start9  end10		11 strand12		13		ID14															overlap15
#1       4910    4927    ACTCAAATAATGAACCGG      59      -       1       ensembl gene    4854    9652    .       -       .       ID=GRMZM2G059865;Name=GRMZM2G059865;biotype=protein_coding      17
#1       12131   12151   TGTAGGTTGGGTATTTTTGAT   2       -       1       MTEC    repeat  10872   15310   .       -       .       class=II;subclass=1;order=TIR;superfamily=[DTC] CACTA;type=Type II Transposons/TIR;name=ZM_CACTA_42     20
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[14]=~/;family=([^;]+)/;
	$family=$1;
	$len=length($tmp[3]);
	$reads_norm=$tmp[4]/$loci{$tmp[3]};
	$reads_norm{$1}+=$reads_norm;
	$reads_len{$1}{21}+=$reads_norm if $len==21;
	$reads_len{$1}{22}+=$reads_norm if $len==22;
	$reads_len{$1}{24}+=$reads_norm if $len==24;
	if ($start!=$tmp[9] && $end!=$tmp[10]) {
		$start=$tmp[9];
		$end=$tmp[10];
		$length=$end-$start+1;
		$length{$family}+=$length;
	}
}
close INPUT;
open OUTPUT,">$opts{o}";
#repeat_family	RPM		RPKM	21-RPM	22-RPM	24-RPM
#Jittery		1.35    0.16    0.00    0.00    1.35
foreach  (sort keys %reads_norm) {
#	if (defined($reads_norm{$_})) {
#		$rpkm{$_}=$reads_norm{$_}*1000000*1000/$total/$len{$_};
#		printf OUTPUT "$_\t%.2f\n",$rpkm{$_};
#	}else{
#		print OUTPUT "$_\t0\n";
#	}
	$rpm{$_}=$reads_norm{$_}*1000000/$total;
	$rpm_len{$_}{21}=$reads_len{$_}{21}*1000000/$total;
	$rpm_len{$_}{22}=$reads_len{$_}{22}*1000000/$total;
	$rpm_len{$_}{24}=$reads_len{$_}{24}*1000000/$total;
	$length{$_}=1000 unless defined($length{$_});
	$rpkm{$_}=$rpm{$_}*1000/$length{$_};
	printf OUTPUT "$_\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",$rpm{$_},$rpkm{$_},$rpm_len{$_}{21},$rpm_len{$_}{22},$rpm_len{$_}{24};
}

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
Usage:    intersect_bedgff_repeat.pl -i gff_cluster -a intersect_gff_result -d gene_description -o result
Function: Compute repeat family region sRNA RPM, also split calculate 21, 22, 24-nt length type
Command:  -i inersect_bed result
          -o cluster_rpkm
          -d cluster
		  -t loci database
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-08-9
Notes:    Originate from intersect_bedgff_repeat_rpkm_length.pl
\n/
    )
}