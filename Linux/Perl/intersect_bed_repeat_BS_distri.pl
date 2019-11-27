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
#chrom0		Start1 End2			name3			BS_%4 strand5 	chr12-6		start13-7 end14-8 geneID15-9	type16-10 strand17-11 overlap18-12 
#3       218912349       218912350       CHH     0.0     G       3       218912213       218914212       ath-rep642857   upstream        +       1
#3       218912350       218912351       CHH     0.0     G       3       218912213       218914212       ath-rep642857   upstream        +       1
#3       218912351       218912352       CHH     0.0     G       3       218912213       218914212       ath-rep642857   upstream        +       1
#3       218912352       218912353       CHH     0.0     G       3       218912213       218914212       ath-rep642857   upstream        +       1
#3       218912370       218912371       CHH     0.0     G       3       218912213       218914212       ath-rep642857   upstream        +       1
#3       218912371       218912372       CHH     0.0     G       3       218912213       218914212       ath-rep642857   upstream        +       1
#3       218912373       218912374       CG      1.0     G       3       218912213       218914212       ath-rep642857   upstream        +       1
	chomp;
	my @tmp=split/\t/;
	$middle=$tmp[2];
	if ($tmp[10] eq 'upstream') {
		if ($tmp[11] eq "+") {
			$pos=2000+1-abs($tmp[8]-$middle);
		}else{
			$pos=2000+1-abs($tmp[7]-$middle);
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos>$i && $pos<=$i+$bin) {
				$count{upstream}{$i}++;
				$BS{upstream}{$i}+=$tmp[4];
				last;
			}
		}
	}elsif($tmp[10] eq 'downstream') {
		if ($tmp[11] eq "+") {
			$pos=abs($tmp[7]-$middle);
		}else{
			$pos=abs($tmp[8]-$middle);
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos>$i && $pos<=$i+$bin) {
				$count{downstream}{$i}++;
				$BS{downstream}{$i}+=$tmp[4];
			last;
			}
		}
	}elsif($tmp[10] eq 'gene') {
		if ($tmp[11] eq "+") {
			$pos=abs($middle-$tmp[7]+1)/($tmp[8]-$tmp[7]+1)*2000; # normalize gene to 2k, relation position
			$pos1=abs($middle-$tmp[7]+1); # position to TSS
			$pos2=2000-$tmp[8]+$middle; # position to TTS
		}else{
			$pos=abs($middle-$tmp[8]+1)/($tmp[8]-$tmp[7]+1)*2000; # normalize gene to 2k, relation position
			$pos1=abs($tmp[8]-$middle+1); # position to TSS
			$pos2=2000+$tmp[7]-$middle; # position to TTS
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos>$i && $pos<=$i+$bin) {
				$count{gene}{$i}++;
				$BS{gene}{$i}+=$tmp[4];
				last;
			}
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos1>$i && $pos1<=$i+$bin) {
				$count{TSS}{$i}++;
				$BS{TSS}{$i}+=$tmp[4];
				last;
			}
		}
		for ($i=$opts{s};$i<$opts{e};$i+=$bin) {
			if ($pos2>$i && $pos2<=$i+$bin) {
				$count{TTS}{$i}++;
				$BS{TTS}{$i}+=$tmp[4];
				last;
			}
		}
	}
}
close INPUT;
open OUTPUT,">$opts{o}";
foreach  $id(@name) {
	for  ($i=$opts{s};$i<$opts{e};$i+=$bin) {
		$count{$id}{$i}=1 unless defined($count{$id}{$i});
		$BS{$id}{$i}=0 unless defined($BS{$id}{$i});
		$BS=$BS{$id}{$i}/$count{$id}{$i};
		printf OUTPUT "$id$i\t%.4f\n",$BS;
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
Usage:    intersect_bed_repeat_BS_distri.pl -d gff_cluster -i intersect_gff_result -t loci_database -o result
Function: use intersect result get assigned region RPM
Command:  -i inersect_bed result
          -o cluster_rpm
          -d cluster
		  -t loci database
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-07-06
Notes:    
\n/
    )
}