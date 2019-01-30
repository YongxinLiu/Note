#! /usr/bin/perl
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;
use Text::NSP::Measures::2D::Fisher::twotailed;



###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'a:b:c:d:', \%opts );
&usage unless ( exists $opts{a} && exists $opts{b} );

$pvalue = &Fisher_Test($opts{a}, $opts{a}+$opts{b}, $opts{a}+$opts{c},$opts{a}+$opts{b}+$opts{c}+$opts{d});
print $pvalue,"\n";

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
Usage:    pvalue.pl -a sample1 -b sample_total -c background1 -d background_total
Function: calculate p-value by fisher test from package Text::NSP::Measures::2D::Fisher::twotailed
Command:  -a sample1
          -b sample_total
          -c background1
		  -d background_total
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-03-19
Notes:    
\n/
    )
}
