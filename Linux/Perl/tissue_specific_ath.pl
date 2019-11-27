#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Tissue number is $opts{d}\n" if defined($opts{d});
$opts{t}=0.5 unless defined($opts{t});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#ID	RPM1	RPM2	RPM3	RPM4	RPM5	RPM6 standard deviation7-12
#ath-SRC9982	1.7	12.45	2.515	12.07	7.25	0.29	1.964051934	6.434671709	1.902117241	0	9.454017664	0
open OUTPUT,">$opts{o}";
#ID	ratio1-6	entropy7	specific8
#ath-SRC9982 0.1	0.1	0.5		0.1	0.1	0.1	0.0001	1
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$total=0;
	my @ratio;
#	print OUTPUT $tmp[0],"\t";
	for $i(1..$opts{d}) { #add 0.01 for exclude 0
		$tmp[$i]+=0.01;
	}
	for $i(1..$opts{d}) { #calculate total RPM
		$total+=$tmp[$i];
	}
	for $i(1..$opts{d}) { #calculate samples ratio
		$ratio[$i]=$tmp[$i]/$total;
#		printf OUTPUT "%.2f\t",$ratio[$_];
	}
	my $sign="";
	for $i(1..$opts{d}) { #calculate samples ratio
		if ($ratio[$i]>0.5 && $tmp[$i]>5 && $tmp[$i+$opts{d}]<$opts{t}*$tmp[$i]){
			$sign="$i";
			last;
		}
#		if ($ratio[$i]>0.25 && $tmp[$i]>5){
#			$sign.="$i\,";
#		}
	}
	$entropy=0;
	my @entropy;
	for $i(1..$opts{d}) { #calculate samples ratio
		$entropy[$i]=$ratio[$i]*log($ratio[$i])/log(2);
		$entropy+=$entropy[$i];
	}
	$entropy*=(-1);
	print OUTPUT "$_\t$entropy\t$sign\n";

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
sub usage {
    die(
        qq/
Usage:    tissue_specific_ath.pl -i inpute_file -o output_file
Function: clacluate shannon entropy and tissue specific SRCs based on tissue average RPM and standard variance 
Command:  -i inpute file name (Must)
          -o output file name (Must)
		  -d tissue number
		  -t samples variation < $opts{t} average, default <0.5
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-10
Notes:    
\n/
    )
}