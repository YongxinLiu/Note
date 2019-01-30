#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Bio::SeqIO;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Main text.
###############################################################################

###############################################################################
#Read the genome in hash
###############################################################################
open DATABASE,"<$opts{d}";
my %database;
my $chr;
while (<DATABASE>) {
	chomp;
	if (/>/) {
		$chr=(split('\s',$_))[0];
		$chr=~s/>//;
#		print $chr,"\n";
	}else{
		$database{$chr}.=$_;
	}
}
close DATABASE;
#print substr($database{$chr},0,100),"\n",substr($database{$chr},-100,100);

open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
#open OUTPUTS,">>$opts{s}";
my %sequence;
while (<INPUT>) {
	chomp;
	my @line=split('\t',$_);
	my $m=scalar(@line);
#Chr1    +       AT1G01010.1     3914    3995    4277    4485    4606    4705    5096    5173    5327    5438
#Chr1    -       AT1G01020.1     8570    8465    8416    8326    8235    7988    7941    7836    7761    7650
	if (defined($line[3])){
#	print $line[3];
	$line[0]=~s/Chr//;
	foreach  (keys %database) {
		if ($_=~m/^$line[0]/i) {
			my $n=3;
			if ($line[1] eq "+"){
				while ($n<$m){
					my $seq1=substr($database{$_},$line[$n]-21,40);
					my $seq2=substr($database{$_},$line[$n+1]-20,40);
					my $site=$line[0].$line[$n];
					if (defined($sequence{$site})){
					}else{
						print OUTPUT ">$line[0]\t$line[1]\t$line[2]\t$line[$n]\n$seq1\n";
						$sequence{$site}=$seq1;
					}
					$site=$line[0].$line[$n+1];
					if (undef($sequence{$site})){
					}else{	
						print OUTPUT ">$line[0]\t$line[1]\t$line[2]\t$line[$n+1]\n$seq2\n";
						$sequence{$site}=$seq2;
					}
					$n=$n+2;
				}
			}else{
				while ($n<$m){
					my $seq1=substr($database{$_},$line[$n]-20,40);
					$seq1= &revcom($seq1);
					my $seq2=substr($database{$_},$line[$n+1]-21,40);
					$seq2= &revcom($seq2);
					my $site=$line[0].$line[$n];
					if (defined($sequence{$site})){
					}else{
						print OUTPUT ">$line[0]\t$line[1]\t$line[2]\t$line[$n]\n$seq1\n";
						$sequence{$site}=$seq1;
					}
					$site=$line[0].$line[$n+1];
					if (defined($sequence{$site})){
					}else{
						print OUTPUT ">$line[0]\t$line[1]\t$line[2]\t$line[$n+1]\n$seq2\n";
						$sequence{$site}=$seq2;
					}	
					$n=$n+2;
				}
					
			}
		}
	}
	}
}
close INPUT;
close OUTPUT;

sub revcom{
	my ($seq)=@_;
	$seq=reverse $seq;
	$seq=~tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

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
Usage:    template.pl -i inpute_file -o output_file
Function:
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-09-20
Notes:    
\n/
    )
}
