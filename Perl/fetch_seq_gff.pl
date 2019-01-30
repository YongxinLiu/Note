#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:l:r:s:t:m:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{l}=0 unless defined($opts{l});
$opts{r}=0 unless defined($opts{r});
$opts{s}=0 unless defined($opts{s});
$opts{m}=1 unless defined($opts{m});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";#reads fasta in hash
#>1
#GAATTCCAAAGCCAAAGATTGCATCAGTTCTGCTGCTATTTCCTCCTATCATTCTTTCTG
#ATGTTGAAAATGATATTAAGCCTAGGATTCGTGAATGGGAGAAGGTATTTTTGTTCATGG
my %database;
my $chr;
while (<DATABASE>) {
	chomp;
	if (/>/) {
		$_=~/>(\w+)/;
		$chr=$1;
#		print $chr,"\n";
	}else{
		$database{$chr}.=$_;
	}
}
close DATABASE;
###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chrom0 source1			type2			Start3	 End4	score5	strand6				description8
#Chr1    omicslab        sRNA_cluster    11017   11176   18.80   +       .       ID=SRC1;length_type=24-nt
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
#	next if $tmp[5]<$opts{s};
#	next if $loci{$tmp[3]}>$opts{m};
	if ($tmp[6] eq '+') {
		$tmp[3]-=$opts{l};
		$tmp[4]+=$opts{r};
	}else{
		$tmp[3]-=$opts{r};
		$tmp[4]+=$opts{l};
	}
	if ($tmp[3]<1) {#treatment start <1, as the 1
		$tmp[3]=1;
	}
	my $seq=substr($database{$tmp[0]},$tmp[3]-1,$tmp[4]-$tmp[3]+1);
	$seq=&revcom($seq) if $tmp[5] eq '-';
	$tmp[8]=~/ID=([\w-]+)\;/;
	print OUTPUT ">$1\n$seq\n";
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
sub revcom{
	my $revcom=$_[0];
	$revcom=reverse $revcom;
	$revcom=~tr/AGCT/TCGA/;
	return $revcom;
}

sub usage {
    die(
        qq/
Usage:    fetch_seq_gff.pl -i file.gff -o file.fasta -d genome.fasta
Function: fetch sequence based on gff file from genome or fasta file
Command:  -i gff file (Must)
          -o output fasta (Must)
          -d fasta database
          -l upstream scale, default=0
          -r downstream scale, default=0
#          -s score scale, default=0
          -t loci data
          -m loci thrshold default=1
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-23
Notes:    Debug problem if start subtract flank less than 1, ouput chromosome end
\n/
    )
}