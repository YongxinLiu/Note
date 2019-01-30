#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:n:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open LOCI,"<$opts{t}";
##sequences0		reads1	loci2
##AAAAAAAAAAAAAACTCGGCC   1       1
#my (%loci,$total);
#while (<LOCI>) {
#	chomp;
#	my @tmp=split/\t/;
#	$loci{$tmp[0]}=$tmp[2];#preserve all the loci number in hash
#	$total+=$tmp[1];
#}
#close LOCI;
open DATABASE,"<$opts{d}";
#AT1G01010
#AT1G01020
#AT1G03987
my @database;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @database,$tmp[0]
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
while (<INPUT>) {
#chrom0 Start1 End2		GI3			Name4 strand5 chrom6	Start7		End8 Ref_count9  Alt_count10 overlap11
#1       3631    5899    AT1G01010       NAC001  +       1       4425    4426    1       0       1
#1       3631    5899    AT1G01010       NAC001  +       1       4647    4648    1       0       1
#1       6788    9130    AT1G01020       ARV1    -       1       8190    8191    3       0       1
	chomp;
	my @tmp=split/\t/;
#	print "$tmp[9]\t$tmp[10]\n";
	$count_ref{$tmp[3]}+=$tmp[9];
	$count_alt{$tmp[3]}+=$tmp[10];
}
close INPUT;
open OUTPUT,">$opts{o}";

printf OUTPUT "GI\t$opts{n}Ref\t$opts{n}Alt\n";

foreach  (@database) {
	if (defined($count_alt{$_})) {
		printf OUTPUT "$_\t$count_ref{$_}\t$count_alt{$_}\n";
	}else{
		print OUTPUT "$_\t0\t0\n";
	}
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
sub usage {
    die(
        qq!
Usage:    intersect_bed_rpm.pl -d GI_list -i intersect_bedVSbed -n name -o output
Function: use intersect result get gene region count
Command:  -i inersect_bed result
          -o gene_count
          -d gene list
		  -n sample name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/5/29
Notes:    
\n!
    )
}