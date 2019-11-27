#!/home/hdsu/software/perl/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
use Text::NSP::Measures::2D::Fisher::twotailed;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:d:o:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{d} );
my $start_time=time;
$opts{t}=1 unless defined($opts{t});
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Sample ID is $opts{i}\ncontrol ID is $opts{d}\n";


###############################################################################
my (%sample,%database);
my @total;
open DATABASE,"<$opts{d}";
<DATABASE>;
	while (<DATABASE>) {
	chomp;
#0gene_id 1transcript_id(s)        2length  3effective_length        4expected_count  5TPM     6FPKM
#TR10000|c0_g1   TR10000|c0_g1_i1,TR10000|c0_g1_i2,TR10000|c0_g1_i3      1290.14 1017.03 16210.00        185.52  446.87
#TR10000|c0_g2   TR10000|c0_g2_i1        291.00  27.94   0.00    0.00    0.00
#TR10000|c0_g3   TR10000|c0_g3_i1        281.00  21.63   0.00    0.00    0.00
#0transcript_id   1gene_id length  2effective_length        3expected_count  4TPM     5FPKM    6IsoPct
#TR10000|c0_g1_i1        TR10000|c0_g1   1362    1088.89 6146.44 65.70   158.26  35.42
#TR10000|c0_g1_i2        TR10000|c0_g1   1357    1083.89 6892.56 74.02   178.29  39.90
#TR10000|c0_g1_i3        TR10000|c0_g1   1079    805.89  3171.00 45.80   110.32  24.69
	my @tmp=split/\t/;
	$database{$tmp[0]}{rpm}=$tmp[6];
	$database{$tmp[0]}{reads}=int($tmp[4]+0.5); # reads must integer
	$total[1]+=$tmp[4];
}
close DATABASE;
open DATABASE,"<$opts{i}";
<DATABASE>;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$sample{$tmp[0]}{rpm}=$tmp[6];
	$sample{$tmp[0]}{reads}=int($tmp[4]+0.5);
	$total[0]+=$tmp[4];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
#ID			rpm1	rpm2	fold	p-value
#SRC1		1		2		0.5		0.001

#Compute P-value by Fisher test
foreach $id(keys %sample) {
	$sample{$id}{fold}=($sample{$id}{rpm}+1)/($database{$id}{rpm}+1);#normalize fold change
#	my $p_value = &Fisher_Test($reads[0], $reads[0]+$total[0], $reads[0]+$reads[1], $reads[0]+$reads[1]+$total[0]+$total[1]);
		# $reads[0]: expression level of geneA in sample1;
		# $total[0]: total mapped reads in sample1;
		# $reads[1]: expression level of geneA in sample2;
		# $total[1]: total mapped reads in sample2;
	$sample{$id}{pvalue} = &Fisher_Test($sample{$id}{reads}, $sample{$id}{reads}+$total[0], $sample{$id}{reads}+$database{$id}{reads}, $sample{$id}{reads}+$database{$id}{reads}+$total[0]+$total[1]);
} 

#Compute FDR
@srcid=sort {$sample{$a}{pvalue}<=>$sample{$b}{pvalue};}  keys %sample; #sort hash by P-value, get the P-value rank
$rank=1;
$length=$#srcid;
foreach $id(@srcid) {
	$sample{$id}{FDR}=$sample{$id}{pvalue}*$length/$rank;
	$rank++;
}

#Output result
foreach my $id(sort keys %sample) {
	printf OUTPUT "$id\t$sample{$id}{rpm}\t$database{$id}{rpm}\t%.2f\t%.12f\t%.8f\n",$sample{$id}{fold},$sample{$id}{pvalue},$sample{$id}{FDR};
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
sub read_rpm{ #reads file rpm in hash
	$ID=$_[0];
	my %list;
	open INPUT,"<$ID";
	while (<INPUT>) {
		chomp;
		@tmp=split/\t/;
		$list{$tmp[0]}=$tmp[1];
	}
	return %list;
	close INPUT;
}

sub compute_pvalue{ #compute p-value
	($n1,$n2,$x1,$y1)=@_;
	print "$x1\t$y1\n";
	return (($n2/$n1)**$y1)*factorial($x1+$y1)/(factorial($x1)*factorial($y1)*(1+$n2/$n1)**($x1+$y1+1));#(N2/N1)**y*(x+y)!/(x!*y!*(1+N2/N1)**(x+y+1))
}

sub factorial{ #compute factorial
        return ($_[0] <= 1) ? 1 : $_[0] * factorial($_[0] - 1);
}

sub Fisher_Test {
	(my $n11, my $n1p, my $np1, my $npp) = @_;
	my $twotailed_value = calculateStatistic( 
						n11=>$n11,
						n1p=>$n1p,
						np1=>$np1,
						npp=>$npp);
	if( (my $errorCode = getErrorCode())){
		print STDERR $errorCode." - ".getErrorMessage();
	}else{
		return $twotailed_value;
	}
}

sub usage {
    die(
        qq/
Usage:    compare_different_trinity.pl -i special -d WT -o output_file
Function: compute trinity genes\/isoforms FPPM fold change, for reducing the low RPM caused big change or divide by zero, fold=treatment+1\/control+1
          calculate p-value by fisher test from package Text::NSP::Measures::2D::Fisher::twotailed
          calculate FDR\/q-value
Command:  -i treatment samples
          -d control, control ID, support two samples average as control
          -o output
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-03-07
Notes:    Update from compare_different_src.pl, only for trinity RSEM.genes.result
\n/
    )
}
