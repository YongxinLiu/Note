#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{s}=2000 unless defined($opts{s});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}"; #selected SRCs list
my %database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
#	print $tmp[0],"\n";
	$database{$tmp[0]}=1;
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my (%overlap,%flank);
while (<INPUT>) {
#ID0			Chr1	Start2	End3	Strand4	RPM5 Strand_ratio6	Copies7	20-24(8-12)						Length_type13	Upstream14 Overlap15 Downstream16
#ath-SRC1        Chr1    1       201     +       9.06    0.60    1.04    0.00    0.07    0.03    0.25    0.65    24-nt   ;       ath-rep1|Chr1|1|107|-|107|Repbase13|DNA|ATREP18;        AT1G01010|Chr1|3631|5899|+|3429|TAIR10|protein_coding_gene|NAC domain containing protein 1;ath-rep2|Chr1|1064|1098|+|862|Repbase13|Simple_repeat|(CACCCCC)n
	chomp;
	my @tmp=split/\t/;
#	print $tmp[0],"\n";
	next unless defined($database{$tmp[0]});#whether this SRCs is our needed
#	print $tmp[15],"\n";
	my @tmp_over=split/\;/,$tmp[15];
	foreach  $over(@tmp_over) {
		if ($over=~/(AT\wG\d+)/) {
			my $id=$1;
			my @tmp1=split/\|/,$over;
			if (defined($overlap{$id})) {
				$tmp1[8]="" unless defined($tmp1[8]);
				$overlap{$id}=$tmp1[8] if length($tmp1[8])>length($overlap{$id});
			}else{
				$overlap{$id}=$tmp1[8];
			}
		}
	}
	&flank_gene($tmp[14]);
	&flank_gene($tmp[16]);
}

foreach  (keys %overlap) {
	print OUTPUT "$_\t$overlap{$_}\n";
}
	print OUTPUT "\n";

foreach  (keys %flank) {
	print OUTPUT "$_\t$flank{$_}\n";
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
sub flank_gene{
	my $flank_anno=$_[0];
	if ($flank_anno=~/(AT\wG\d+)/) {
		my @tmp=split('\|',$flank_anno);
		if ($tmp[5]<$opts{s}) {
			my @tmp2=split/\;/,$tmp[8];
			$flank{$tmp[0]}=$tmp2[0];
		}
	}

}
sub usage {
    die(
        qq/
Usage:    SRCs_geneID_list_ath.pl.pl -i inpute_file -o output_file
Function: Get SRCs overlapped gene ID list for Gene Ontology (GO) analysis
Command:  -i SRCs annotation (Must)
          -o Gene list (Must)
          -d SRC list (Must)
		  -s flank scale, default<2000bp
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-04
Notes:    
\n/
    )
}