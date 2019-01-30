#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:n:', \%opts );
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

$total=$#database-$opts{n};



###############################################################################
#Main text.
###############################################################################
my %hash;
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my $count=0;
my $column;
while (<INPUT>) {
	my @sum;
	$count++;
	chomp;
	$line=int($_*$total);
	while (defined($hash{$line})) {
#		print $line,"\n";
		$line=int($line/2+1);
	}
	$hash{$line}=0;
	for ($i=$line;$i<$line+$opts{n} ;$i++) {
		chomp($database[$line]);
		my @tmp=split/\t/,$database[$line];
#		print  "$tmp[1]\n";
		$column=$#tmp;
		for (1..$#tmp) {
			$sum[$_]+=$tmp[$_];
		}
	}
	print OUTPUT "Rand$count";
#	print $column,"\n";
	for (1..$column) {
#		print $sum[$_];
		$sum[$_]=$sum[$_]/$opts{n};
		print OUTPUT "\t$sum[$_]";
	}
	print OUTPUT "\n";
	
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
Usage:    random_line_seed1.pl -i random_seed -o output_file -d database
Function: base on random seed, select random N lines from database
Command:  -i randome_seed text (Must)
          -o output file name
          -d database file name
          -h database header line number, default 0
          -n number
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-07
Notes:    
\n/
    )
}