#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{d} );
$opts{d}="result_k1-c/database.txt" unless defined($opts{d});
#$opts{o}="$opts{i}.xls" unless defined($opts{o});
#`rm $opts{i}.xls`;
$opts{o}="$opts{i}.xls";
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my %database; #database in hash
$header=<DATABASE>;
#print $header;
while (<DATABASE>) {
#	print $_;
	my @tmp=split/\t/;
#	print $tmp[0],"\n";
	$database{$tmp[0]}=$_;
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
print OUTPUT "\t$header";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT "\n$tmp[0]\n$tmp[1]\n";
	if ($tmp[1]>0) {
		@otu=split(/,/,$tmp[2]);
		foreach (@otu) {
			print  $database{$_},"\n";
			print OUTPUT $database{$_};
		}
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
Usage:    vennNumAnno.pl -i inpute_file -o output_file -d database	-h header num
Function: Annotate vennNumberGenerator result
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017-7-8
Notes:    
\n/
    )
}
