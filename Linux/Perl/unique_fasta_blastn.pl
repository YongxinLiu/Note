#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:c:s:', \%opts );
&usage unless ( exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{c}=97 unless defined($opts{c});
$opts{s}=99 unless defined($opts{s});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{i}";#reads fasta in hash
#>1
#GAATTCCAAAGCCAAAGATTGCATCAGTTCTGCTGCTATTTCCTCCTATCATTCTTTCTG
#ATGTTGAAAATGATATTAAGCCTAGGATTCGTGAATGGGAGAAGGTATTTTTGTTCATGG
my %database;
my %len;
my $chr;
while (<DATABASE>) {
	chomp;
	if (/>/) {
		$_=~/>([\w\.\|]+)/;
		$chr=$1;
#		print $chr,"\n";
	}else{
		$database{$chr}.=$_;
		$len{$chr}=length($database{$chr});
	}
#	print "$chr\t$len{$chr}\t$database{$chr}\n";
}
close DATABASE;
my %unique;
###############################################################################
#Main text.
###############################################################################
#chrom0 source1			type2			Start3	 End4	score5	strand6				description8
#Chr1    omicslab        sRNA_cluster    11017   11176   18.80   +       .       ID=SRC1;length_type=24-nt
open INPUT,"<$opts{d}";
open OUTPUT,">$opts{o}";
my $count=0;
my $query="F";
my %seq;
my %id;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
#qseqid0            sseqid1	pident2 length3 mismatch4 gapopen5 qstart6 qend7    sstart8     send9           evalue10 bitscore11
#IGS-A-long      scaffold51788   92.67   2700    173     23      122     2810    7463    10148   0.0     3866
#IGS-A-long      scaffold51788   91.41   1641    123     18      122     1756    17186   18814   0.0     2233
	next if $tmp[2]<$opts{s}; # remove similarity < 99%
#	print "$len{$tmp[0]}\t";
	next if $tmp[3]/$len{$tmp[0]}*100<$opts{c}; # remove coverage < 97%
	next if !defined($database{$tmp[0]}); # remove used reads
#	print $tmp[3]/$len{$tmp[0]}*100,"\t";
#	print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\n";
	if ($query ne $tmp[0]) {
		delete($database{$query});
		$query=$tmp[0];
		$count++;
		$id{$count}=$tmp[0];
	}
	if ($len{$tmp[0]}>$len{$tmp[1]}) {
		$seq{$count}=$database{$tmp[0]};
	}else{
		$seq{$count}=$database{$tmp[1]};
	}
	next if $tmp[0] eq $tmp[1]; # remove self
	$id{$count}.="_".$tmp[1];
	delete($database{$tmp[1]});
}

foreach  (1..$count) {
	print OUTPUT ">$id{$_}\n$seq{$_}\n";
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
sub revcom{
	my $revcom=$_[0];
	$revcom=reverse $revcom;
	$revcom=~tr/AGCT/TCGA/;
	return $revcom;
}

sub usage {
    die(
        qq/
Usage:    unique_fasta_blastn.pl -i raw.fasta -o file.fasta -d raw_self.blastn
Function: remove redandency seq from sequencing result by coverage and similarity
Command:  -i raw fasta (Must)
          -o output fasta (Must)
          -d self blastn
          -c coverage, default > 97%
          -s similarity, default > 99%
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-10-26
Notes:    
\n/
    )
}
