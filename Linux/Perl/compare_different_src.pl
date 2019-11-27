#!/usr/bin/perl -w
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
# ID	count	RPM
#57923_655	135	5.82444099388261
#29173_595	2	0.0862880147241868
#29318_306	3	0.12943202208628

while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}{rpm}=$tmp[1];
	$database{$tmp[0]}{reads}=$tmp[2];
	$total[1]+=$tmp[2];
}
close DATABASE;
open DATABASE,"<$opts{i}";
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$sample{$tmp[0]}{rpm}=$tmp[1];
	$sample{$tmp[0]}{reads}=$tmp[2];
	$total[0]+=$tmp[2];
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
Usage:    compare_different_src.pl -i sample1 -d sample2 -o output file
Function: compute SRCs RPM fold change, for reducing the low RPM caused big change or divide by zero, fold=treatment+1\/control+1
          calculate p-value by fisher test from package Text::NSP::Measures::2D::Fisher::twotailed
          calculate FDR\/q-value
Command:  -i treatment samples
          -d control, control ID, support two samples average as control
          -o output
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-30
Notes:    
\n/
    )
}
