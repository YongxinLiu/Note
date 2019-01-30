#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}"; # read copies number
my %database; #database in hash
while (<DATABASE>) {
#chr start end	copies
#1	2	3	2
#1	3	4	2
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}{$tmp[1]}=$tmp[3];
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#>1_0_23
#GAATTCCAAAGCCAAAGATTGCA
	chomp;
	if (/>/) {
		$_=~/>([\w_]+)_(\d+)_(\d+)/;
		$chr=$1; # 内部变量需保存，在子过程中失效
		$pos=$2;
		$seq=<INPUT>;
		chomp($seq);
		next unless ($seq=~/^CC/ || $seq=~/GG$/);
		next unless $database{$chr}{$pos}<2;
		#print $database{$chr}{$pos},"\n";
		print OUTPUT "$chr\t$pos\t",$pos+23,"\t$seq\t$database{$chr}{$pos}\n";
	}
	
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
Usage:    CRISPR_target_copy.pl -i inpute_file -o output_file -d database	-h header num
Function: get CRISPR target, and copies < 4
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}