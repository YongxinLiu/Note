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
my %reads;
open DATABASE,"<$opts{d}";
#chrom0			Start1 End2			name3		score4 strand5		Chr6		start7	 end8	ID9	overlap10
#Bd1	37013338	37013491	CentBd	6	-	Bd1	37013330	37013351	4.38503	13
#Bd1	37013338	37013491	CentBd	6	-	Bd1	37013351	37013359	3.65419	8
#Bd1	37013338	37013491	CentBd	6	-	Bd1	37013359	37013367	4.38503	8
#Bd1	37013338	37013491	CentBd	6	-	Bd1	37013367	37013369	5.11586	2
my %score;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$id="$tmp[0]\t$tmp[1]\t$tmp[2]\t";
	if (defined($score{$id})) {
		$score{$id}=$tmp[10] if $tmp[10]>$score{$id};
	}else{
		$score{$id}=$tmp[10];
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chrom0			Start1 End2			name3		score4 strand5
#Bd1	37013338	37013491	CentBd	6	-
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$id="$tmp[0]\t$tmp[1]\t$tmp[2]\t";
	if (defined($score{$id})) {
		print OUTPUT "$id",$tmp[2]-$tmp[1],"\t$score{$id}\t$tmp[5]\n";
	}else{
		print OUTPUT "$id",$tmp[2]-$tmp[1],"\t0\t$tmp[5]\n";
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
Usage:    intersect_bed_count.pl -i bed_list -d intersectbed -o output -t total reads for RPM
Function: get the intersect value max socre for each line
Command:  -i inpute file name (Must)
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
