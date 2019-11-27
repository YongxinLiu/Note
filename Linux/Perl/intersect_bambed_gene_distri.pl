#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:s:e:b:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{s}=0 unless defined($opts{h});
$opts{e}=2000 unless defined($opts{h});
$opts{b}=50 unless defined($opts{h});
$bin=$opts{e}/$opts{b};

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
my @name=("upstream","TSS","TTS","downstream","upstream","gene","downstream");
while (<INPUT>) {
	my ($middle,$pos);
#chrom0		Start1 End2			name3					mQ4 strand5 6	7	8		9	10	11	chr12		start13 end14 geneID15	type16 strand17 overlap18 
#C159934294	4	102	HWI-ST958:140:8:2108:7953:127341#0/1	6	+	4	102	0,0,0	1	98,	0,	C159934294	1	115 EPlTURG00000002234	upstream	+	98
#C159934294^I4^I102^IHWI-ST958:140:8:2108:7953:127341#0/1^I6^I+^I4^I102^I0,0,0^I1^I98,^I0,^IC159934294^I1^I115^IEPlTURG00000002234^Iupstream^I+^I98$
#C159934294^I7^I89^IHWI-ST833:185:7:2204:1280:30927#0/2^I0^I-^I7^I89^I0,0,0^I1^I82,^I0,^IC159934294^I1^I115^IEPlTURG00000002234^Iupstream^I+^I82$
	chomp;
	my @tmp=split/\t/;
	$middle=($tmp[1]+$tmp[2])/2;
	if ($tmp[16] eq 'upstream') {
		if ($tmp[17] eq "+") {
			$pos=2000+1-abs($tmp[14]-$middle);
		}else{
			$pos=2000+1-abs($tmp[13]-$middle);
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos>$i && $pos<=$i+$bin) {
				$count{upstream}{$i}++;
				last;
			}
		}
	}elsif($tmp[16] eq 'downstream') {
		if ($tmp[17] eq "+") {
			$pos=abs($tmp[13]-$middle);
		}else{
			$pos=abs($tmp[14]-$middle);
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos>$i && $pos<=$i+$bin) {
				$count{downstream}{$i}++;
				last;
			}
		}
	}elsif($tmp[16] eq 'gene') {
		if ($tmp[17] eq "+") {
			$pos=abs($middle-$tmp[13]+1)/($tmp[14]-$tmp[13]+1)*2000; # normalize gene to 2k, relation position
			$pos1=abs($middle-$tmp[13]+1); # position to TSS
			$pos2=2000-$tmp[14]+$middle; # position to TTS
		}else{
			$pos=abs($middle-$tmp[14]+1)/($tmp[14]-$tmp[13]+1)*2000; # normalize gene to 2k, relation position
			$pos1=abs($tmp[14]-$middle+1); # position to TSS
			$pos2=2000+$tmp[13]-$middle; # position to TTS
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos>$i && $pos<=$i+$bin) {
				$count{gene}{$i}++;
				last;
			}
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos1>$i && $pos1<=$i+$bin) {
				$count{TSS}{$i}++;
				last;
			}
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos2>$i && $pos2<=$i+$bin) {
				$count{TTS}{$i}++;
				last;
			}
		}
	}
}
close INPUT;
open OUTPUT,">$opts{o}";
foreach  $id(@name) {
	for  ($i=$opts{s};$i<$opts{e};$i+=$bin) {
		$count{$id}{$i}=0 unless defined($count{$id}{$i});
		print OUTPUT "$id$i\t$count{$id}{$i}\n";
	}
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
Usage:    intersect_bambed_gene_distri.pl -d gff_cluster -i intersect_gff_result -t loci_database -o result
Function: use intersect result get assigned region RPM
Command:  -i inersect_bed result
          -o cluster_rpm
          -d cluster
		  -t loci database
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-07-16
Notes:    
\n/
    )
}