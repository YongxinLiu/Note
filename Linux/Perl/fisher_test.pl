#!/home/yongxin/data38/software/perl_bin/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
use Text::NSP::Measures::2D::Fisher::twotailed;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'a:b:c:d:', \%opts );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

#	my $p_value = &Fisher_Test($reads[0], $reads[0]+$total[0], $reads[0]+$reads[1], $reads[0]+$reads[1]+$total[0]+$total[1]);
		# $reads[0]: expression level of geneA in sample1;
		# $total[0]: total mapped reads in sample1;
		# $reads[1]: expression level of geneA in sample2;
		# $total[1]: total mapped reads in sample2;
$pvalue = &Fisher_Test($opts{a}, $opts{a}+$opts{b}, $opts{a}+$opts{c}, $opts{a}+$opts{b}+$opts{c}+$opts{d});
 
print "P-value = ",$pvalue,"\n";

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
Usage:    fisher_test.pl -a sample1 -b sample1_total -c sample2 -d samples2_total
Function: compute SRCs RPM fold change, for reducing the low RPM caused big change or divide by zero, fold=treatment+1\/control+1
          calculate p-value by fisher test from package Text::NSP::Measures::2D::Fisher::twotailed
Command:  -a sample1
          -b sample1_total
          -c sample2
          -d sample2_total
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-11-25
Notes:    
\n/
    )
}
