#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:g:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{g}=10 unless defined($opts{g});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my %expr; #expression hash
my $i; #count line number
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$expr{$tmp[0]}=$tmp[1];
	$i++;
}

$group=int($i/$opts{g});
print "Total $i line.\tGroup into $opts{g} groups.\tEach group number is $group.\n";

my $groupid=0;
my $count=0;
foreach  (sort {$expr{$a}<=>$expr{$b};}  keys %expr ) {
	$count++; 
	if ($count>$group) { #When count out one group, increment groupID and zero setting count
		$count=1;
		$groupid++;
	}
	$groupid=$opts{g}-1 if ($groupid>=$opts{g}); #all residues into last group
	print OUTPUT "$groupid\t$_\t$expr{$_}\n";
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
Usage:    group_by_expression.pl -i inpute_file -o output_file -g group_number
Function: group gene or sRNA by expression, default from low to high into 10 groups
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
          -g group number, default is 10
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-08-22
Notes:    
\n/
    )
}