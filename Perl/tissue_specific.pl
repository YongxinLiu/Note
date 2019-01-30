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
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{t}=0.5 unless defined($opts{t});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#ID	RPM1	RPM2	RPM3
#SRC1 3.5	5.3	10.8
open OUTPUT,">$opts{o}";
#ID	ratio1	ratio2	ratio3	entropy
#SRC1 0.1	0.1	0.5		0.0001
$head=<INPUT>;
chomp($head);
@tissue=split/\t/,$head;
print OUTPUT "$head\tentropy\ttissue\n";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$total=0;
	my @ratio;
#	print OUTPUT $tmp[0],"\t";
	for $i(1..$#tmp) { #add 0.01 for exclude 0
		$tmp[$i]+=0.01;
	}
	for $i(1..$#tmp) { #calculate total RPM
		$total+=$tmp[$i];
	}
	for $i(1..$#tmp) { #calculate samples ratio
		$ratio[$i]=$tmp[$i]/$total;
#		printf OUTPUT "%.2f\t",$ratio[$_];
	}
	my $sign="";
	for $i(1..$#tmp) { #calculate samples ratio
		if ($ratio[$i]>$opts{t}){
			$sign=$i.$tissue[$i];
			last;
		}
#		if ($ratio[$i]>0.25 && $tmp[$i]>5){
#			$sign.="$i\,";
#		}
	}
	$entropy=0;
	my @entropy;
	for $i(1..$#tmp) { #calculate samples ratio
		$entropy[$i]=$ratio[$i]*log($ratio[$i])/log(2);
		$entropy+=$entropy[$i];
	}
	$entropy*=(-1);
	printf OUTPUT "$_\t%.4f\t$sign\n",$entropy if defined($sign); 
	print OUTPUT "$_\t%.4f\t\n",$entropy unless defined($sign); 

	
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
Usage:    tissue_specific.pl -i inpute_file -o output_file
Function: clacluate tissue specific or universal expressed SRCs, calculate shannon entropy, and select ratio > 0.5 as specific
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -t ratio specific, default 0.5
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-15
Notes:    Includ analysis head line
\n/
    )
}