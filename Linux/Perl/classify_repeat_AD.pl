#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:f:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#ID	AA	DD	AD	AA+DD
#1_9033	167.50	130.97	148.37	149.25	TE	LTR/Copia|ANGELA6_TM-LTR
#2_8779	31.75	102.03	51.64	66.86	Unknown	
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my $tag;
	my @tmp=split/\t/;
	if ($tmp[1]>$opts{f}*$tmp[2]) {
		$tag="AA_Specific";
	}elsif ($opts{f}*$tmp[1]<$tmp[2]){
		$tag="DD_Specific";
	}else{
		$tag="AD_Common";
	}
	print OUTPUT "$tag\t$_\n";
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
Usage:    classify_repeat_AD.pl -i us.value.anno -o us.value.anno.AD
Function: classify AA, DD specific repeat by RPM, AA>4*DD is specifc, other is common
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -f fold change
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-11-05
Notes:    
\n/
    )
}