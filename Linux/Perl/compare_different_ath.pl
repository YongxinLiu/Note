#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:d:o:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{d} );
my $start_time=time;
$opts{t}=1 unless defined($opts{t});
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Numerator ID is $opts{i}\ncontrol ID is $opts{d}\n";


###############################################################################
#Read the database in memory(opt)
my (%control,%control1,%control2);
%sample=&read_rpkm($opts{i});
($c1,$c2)=split(/\+/,$opts{d});#control may be one or two samples
%control1=&read_rpkm($c1);
if (defined($c2)) {
	%control2=&read_rpkm($c2);
	foreach (keys %control1) {
		$control{$_}=($control1{$_}+$control2{$_})/2;
	}
}else{
	%control=%control1;
}

###############################################################################
#Main text.
###############################################################################
$n=(split/\//,$opts{i})[1];
$d=(split/\//,$opts{d})[1];
$output="$opts{o}"."$n"."_"."$d";
print "Output file is $output\n";
open OUTPUT,">$output";
#ID			RPKM	RPKM	fold
#SRC1		1		2		0.5
my $i=0;
foreach my $id(sort keys %sample) {
	if ($control{$id}<$opts{t}) {
		if ($sample{$id}>=$opts{t}) {
			$fold=$sample{$id}/$opts{t};
		}else{
			$fold=1;
		}
	}else{
		$fold=$sample{$id}/$control{$id};
	}
	printf OUTPUT "$id\t$sample{$id}\t$control{$id}\t%.2f\n",$fold;
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
sub read_rpkm{
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



sub usage {
    die(
        qq/
Usage:    compare_different_ath.pl -n numerator -d control 
Function: compute repeat region expression RPKM fold change
Command:  -n numerator, need compute fold change samples ID
          -d control, control ID, support two samples average as control
          -t threshold of RPKM, default=1, number lower can filter as not change
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-11-12
Notes:    
\n/
    )
}