#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:t:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
my %anno;
open DATABASE,"<$opts{d}";
#chrom0			Start1 End2		name3		Chr4		start5	 end6		ID7			type8	strand9	overlap10
#C167343617      1447    2460    AA_18   C167343617      393     4195    TRIUR3_30907    gene    +       1013
#C167381557      5097    5803    AA_22   C167381557      4926    7433    TRIUR3_32192    gene    +       706
#scaffold10001   102619  103500  AA_207  scaffold10001   101337  103336  TRIUR3_11411    upstream        +       717
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($anno{$tmp[3]}{overlap})) {
		next if $anno{$tmp[3]}{overlap}>$tmp[10];
		$anno{$tmp[3]}{id}=$tmp[7];
		$anno{$tmp[3]}{type}=$tmp[8];
		$anno{$tmp[3]}{overlap}=$tmp[10];
	}else{
		$anno{$tmp[3]}{id}=$tmp[7];
		$anno{$tmp[3]}{type}=$tmp[8];
		$anno{$tmp[3]}{overlap}=$tmp[10];
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#scaffold126421  6036    6123    IGS8
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($anno{$tmp[3]}{type})) {
		print OUTPUT "$_\t$anno{$tmp[3]}{type}\t$anno{$tmp[3]}{id}\t$anno{$tmp[3]}{overlap}\n";
	}else{
		print OUTPUT "$_\tintergenic\n";
	}
}
close OUTPUT;
close INPUT;

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
Usage:    intersect_bed_peak_gene_flank.pl -i bed_list -d intersectbed -o output -t total reads for RPM
Function: get the intersect result into peak bed add overlapped type, geneID and overlapped length
Command:  -i inpute peak bed (Must)
          -o output file name (Must)
          -d intersectBed 
          -t total reads
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0 count reads and RPM in each region
Update:   2014-06-10
Notes:    
\n/
    )
}
