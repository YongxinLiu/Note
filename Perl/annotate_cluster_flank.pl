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

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my (%repeat,%gene);
while (<DATABASE>) {
	chomp;
	#Chr1    RepeatMasker    repeat  78393   78473   .       +       .       ID=31;Note=Low_complexity;Annotation=T-rich
	#Chr1    RepeatMasker    repeat  78494   78577   .       +       .       ID=32;Note=Satellite;Annotation=ATCLUST1
	#Chr1    miRbase19       miRNA   78932   78952   .       -       .       ID=ath-miR165a;Note=miRNA;Annotation=ath-miR165a
	my @tmp=split/\t/;
	$tmp[8]=~/ID=([^;]+);Note=([^;]+);Annotation=([^;]*)/; #some anno is blank, need match 0+, use * instead +
	($id,$anno_type,$description)=($1,$2,$3);
#	print "$1\t$2\t$3\t$tmp[8]\n";
	if ($tmp[2]=~/repeat/) {
#	print "$1\t$2\t$3\t$tmp[0]\n";
		$repeat{$tmp[0]}{$id}{source}=$tmp[1];
		$repeat{$tmp[0]}{$id}{start}=$tmp[3];
		$repeat{$tmp[0]}{$id}{end}=$tmp[4];
		$repeat{$tmp[0]}{$id}{strand}=$tmp[6];
		$repeat{$tmp[0]}{$id}{anno_type}=$anno_type; #use type, types, note clash with reserve word
		$repeat{$tmp[0]}{$id}{description}=$description;
	}else{
		$gene{$tmp[0]}{$id}{source}=$tmp[1];
		$gene{$tmp[0]}{$id}{start}=$tmp[3];
		$gene{$tmp[0]}{$id}{end}=$tmp[4];
		$gene{$tmp[0]}{$id}{strand}=$tmp[6];
		$gene{$tmp[0]}{$id}{anno_type}=$anno_type;
		$gene{$tmp[0]}{$id}{description}=$description;
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	#Chr1    omicslab        sRNA_cluster    8       1025    57.05   -       .       ID=SRC1;length_type=24-nt
	chomp;
	my ($leaf_gene,$leaf_repeat,$right_gene,$right_repeat);#set all the description variation
	my ($dis_lg,$dis_lr,$dis_rg,$dis_rr)=(1000000,1000000,1000000,1000000);#Initial distance 1M
	my @tmp=split('\t',$_);
	my @key2=keys %{$gene{$tmp[0]}};#reture gene ID in same chromosome
	foreach $key (@key2) {
		if ($gene{$tmp[0]}{$key}{end}<$tmp[3]) {#find the leaf nearest gene
			$distance=$tmp[3]-$gene{$tmp[0]}{$key}{end}-1;
			if ($distance<$dis_lg) {
				$dis_lg=$distance;
				$leaf_gene="$key\|$tmp[0]\|$gene{$tmp[0]}{$key}{start}\|$gene{$tmp[0]}{$key}{end}\|$gene{$tmp[0]}{$key}{strand}\|$distance\|$gene{$tmp[0]}{$key}{source}\|$gene{$tmp[0]}{$key}{anno_type}\|$gene{$tmp[0]}{$key}{description}";
				#ID,chr,start,end,strand,distance,source,type,description
			}
			next;
		}elsif ($gene{$tmp[0]}{$key}{start}>$tmp[4]){#find the right nearest gene
			$distance=$gene{$tmp[0]}{$key}{start}-$tmp[4]-1;
			if ($distance<$dis_rg) {
				$dis_rg=$distance;
				$right_gene="$key\|$tmp[0]\|$gene{$tmp[0]}{$key}{start}\|$gene{$tmp[0]}{$key}{end}\|$gene{$tmp[0]}{$key}{strand}\|$distance\|$gene{$tmp[0]}{$key}{source}\|$gene{$tmp[0]}{$key}{anno_type}\|$gene{$tmp[0]}{$key}{description}";#ID,chr,start,end,strand,distance,source,type,description
			}
			next;
		}
	}

	@key2=keys %{$repeat{$tmp[0]}};#reture repeat ID in same chromosome
	foreach $key (@key2) {
		if ($repeat{$tmp[0]}{$key}{end}<$tmp[3]) {#find the leaf nearest repeat
			$distance=$tmp[3]-$repeat{$tmp[0]}{$key}{end}-1;
			if ($distance<$dis_lr) {
				$dis_lr=$distance;
				$leaf_repeat="$key\|$tmp[0]\|$repeat{$tmp[0]}{$key}{start}\|$repeat{$tmp[0]}{$key}{end}\|$repeat{$tmp[0]}{$key}{strand}\|$distance\|$repeat{$tmp[0]}{$key}{source}\|$repeat{$tmp[0]}{$key}{anno_type}\|$repeat{$tmp[0]}{$key}{description}";#ID,chr,start,end,strand,distance,source,type,description
			}
			next;
		}elsif ($repeat{$tmp[0]}{$key}{start}>$tmp[4]){#find the right nearest repeat
			$distance=$repeat{$tmp[0]}{$key}{start}-$tmp[4]-1;
			if ($distance<$dis_rr) {
				$dis_rr=$distance;
				$right_repeat="$key\|$tmp[0]\|$repeat{$tmp[0]}{$key}{start}\|$repeat{$tmp[0]}{$key}{end}\|$repeat{$tmp[0]}{$key}{strand}\|$distance\|$repeat{$tmp[0]}{$key}{source}\|$repeat{$tmp[0]}{$key}{anno_type}\|$repeat{$tmp[0]}{$key}{description}";#ID,chr,start,end,strand,distance,source,type,description
			}
			next;
		}
	}

	print OUTPUT "$_\t$leaf_gene;$leaf_repeat\t$right_gene;$right_repeat\n";
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
 usage {
    die(
        qq/
Usage:    annotate_cluster_flank.pl -i inpute_file -o output_file -d annotation
Function: annotate cluster left and right  nearest gene and repeat
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-04-05
Notes:    
\n/
    )
}