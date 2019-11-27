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
open DATABASE,"<$opts{d}";
my %database; #database in hash
# EPlTURG00000000001	-	-	EPlTURG00000000001	MIR1122	-	scaffold9027:9977-10070	-	- OK	0	0	0	OK	0	0	0	OK
#1,10,14,18
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}="$tmp[9]\t$tmp[13]\t$tmp[17]\n";
}

close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
# scaffold125401	59	183	EPlTURG00000000616	gene	-
# ID	reads	RPM
# EPlTURG00000002234	0	0
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($database{$tmp[0]})) {
		print OUTPUT "$tmp[0]\t$database{$tmp[0]}";
	}else{
		print OUTPUT "$tmp[0]\t0\t0\t0\n";
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
Usage:    cuffdiff_gene_value.pl -i gene.bed -o gene_list -d cuffdiff	-h header num
Function: get cuffdiff gene expression, and show in AA\/AD\/DA, output order as input
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