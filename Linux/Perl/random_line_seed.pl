#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{h}=0 unless defined($opts{h});

###############################################################################
#Read the data in array
###############################################################################
my @database;
open DATABASE,"<$opts{d}";
while ($opts{h}>0) { #filter header
	<DATABASE>;
	$opts{h}--;
}
while (<DATABASE>) {
	push @database,$_;
}
close DATABASE;

$total=$#database;



###############################################################################
#Main text.
###############################################################################
my %hash;
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	$line=int($_*$total);
	while (defined($hash{$line})) {
		print $line,"\n";
		$line=int($line/2+1);
	}
	$hash{$line}=0;
	print OUTPUT $database[$line];
	
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
Usage:    random_line_seed.pl -i random_seed -o output_file -d database
Function: base on random seed, select random line from database
Command:  -i randome_seed text (Must)
          -o output file name
          -d database file name
          -h database header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-07
Notes:    
\n/
    )
}