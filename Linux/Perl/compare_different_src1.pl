#!/home/yongxin/data/software/perl_bin/bin/perl -w
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
my $i=0;
foreach my $id(sort keys %sample) {
	$sample{$id}{fold}=($sample{$id}{rpm}+1)/($database{$id}{rpm}+1);#normalize fold change
#	my $p_value = &Fisher_Test($reads[0], $reads[0]+$total[0], $reads[0]+$reads[1], $reads[0]+$reads[1]+$total[0]+$total[1]);
		# $reads[0]: expression level of geneA in sample1;
		# $total[0]: total mapped reads in sample1;
		# $reads[1]: expression level of geneA in sample2;
		# $total[1]: total mapped reads in sample2;
	$sample{$id}{pvalue} = &Fisher_Test($sample{$id}{reads}, $sample{$id}{reads}+$total[0], $sample{$id}{reads}+$database{$id}{reads}, $sample{$id}{reads}+$database{$id}{reads}+$total[0]+$total[1]);
} 
foreach my $id(sort keys %sample) {
	printf OUTPUT "$id\t$sample{$id}{rpm}\t$database{$id}{rpm}\t%.2f\t%.7f\n",$sample{$id}{fold},$sample{$id}{pvalue};
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
Usage:    compare_different_src1.pl -n numerator -d control -o output directory
Function: compute SRCs RPM fold change, for reducing the low RPM caused big change or divide by zero, fold=treatment+1\/control+1
          calculate p-value by fisher test from package Text::NSP::Measures::2D::Fisher::twotailed
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
