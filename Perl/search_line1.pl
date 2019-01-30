#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:c:r:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{c}=0 unless defined($opts{c});
$opts{r}=0 unless defined($opts{r});

my %database;
open DATABASE,"<$opts{d}";
while (<DATABASE>) {
	chomp;
	my @tmp=split/\s+/;
	#print $tmp[0],"\n";
	$database{$tmp[$opts{c}]}=0;
#	print "$tmp[$opts{c}]\n";
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
my $i=0;
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	my @tmp=split/\t/;
	if (defined($database{$tmp[$opts{r}]})) {
#		print $database{$tmp[$opts{r}]},"\n";
		print OUTPUT $_;
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
Usage:    search_line1.pl -i line_file -d gene ID -o line_file
Function: search sequence or gene ID from a dataset
Command:  -i line_file (Must)
          -d GI (Must)
          -o searching result
          -c ID column, default=0
          -r ID in database column, default=0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-03-25
Notes:    Add column select GI column in database, add database ID column, using hashing must same ID
\n/
    )
}
