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
#database in hash£¬ read right reads in database
my %database;
while (<DATABASE>) {
	my @tmp=split/\#/;
	$database{$tmp[0]}=$_.<DATABASE>.<DATABASE>.<DATABASE>;
#	print "$tmp[0]\n$database{$tmp[0]}\n";
}

close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUTA,">$opts{o}\_1\.fastq";
open OUTPUTB,">$opts{o}\_2\.fastq";
#open OUTPUTC,">$opts{o}\_uP\.fastq";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) { #read leaf reads by line
	my @tmp=split/\#/;
	$reads=$_.<INPUT>.<INPUT>.<INPUT>;
	if (defined($database{$tmp[0]})) {
		print OUTPUTA $reads;
		print OUTPUTB $database{$tmp[0]};
		delete $database{$tmp[0]};
#	}else{
#		print OUTPUTC $reads;
	}
}

close INPUT;
close OUTPUTA;
close OUTPUTB;
#close OUTPUTC;

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
Usage:    fastq_pair.pl -i leaf_reads.fastq -d right_reads.fastq -o output_file_prefix(leaf, right, unpair)
Function: read filtered pair fastq, output paired fastq and unpaired fastq. First read right in memory as database, then read each line of left; compared with right and 
Command:  -i leaf_reads.fastq (Must)
          -o output_file_prefix (Must)
          -d right_reads.fastq
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-08-22
Notes:    
\n/
    )
}