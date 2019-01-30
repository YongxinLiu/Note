#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i});
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#database in array
my @tmp1;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @tmp1,$tmp[0];
#	print $tmp[0],"\n";
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
my %data;
while (<INPUT>) {
	chomp;
	my @tmp=split/\s+/;
	$data{$tmp[2]}=$tmp[1];
#	print "$tmp[2]\t$tmp[1]\n";
}

close INPUT;
foreach  $a(@tmp1) {
	if (defined($data{$a})) {
#		print $a,"\n";
		print $data{$a},"\t\t",$a,"\n";
	}else{
		print "0\t\t$a\n";
	}
}


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
Usage:    sort_by_dict2.pl -i inpute_file -d dict
Function: sort by given list, and output list value on screen
Command:  -i inpute file name (Must), uniqc -c output
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2012-05-04
Notes:    
\n/
    )
}