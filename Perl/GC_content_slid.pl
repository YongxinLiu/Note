#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:b:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{b}=10000 unless defined($opts{b});


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
#ID	length	GCnumber	GCcontent
my %database;
my %total=(A,0,G,0,C,0,T,0,N,0);
while (<INPUT>) {
	chomp;
	if (/>/) {
		$_=~/>([\w-]+)/;
		$chr=$1;
#		print $chr,"\n";
	}else{
		$database{$chr}.=$_;
	}
}
foreach  (keys %database) {
	my @nt=split//,$database{$_};
	foreach  (@nt) {
		$i++;
		$eachgc{$_}++;
		if ($i==$opts{b}) {
			$j++;
			$gc=$eachgc{G}+$eachgc{C};
			$gc_content=$gc/$opts{b};
			print OUTPUT $opts{b}*$j,"\t$gc_content\n";
			$i=0;
			$eachgc{G}=0;
			$eachgc{C}=0;
		}
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
Usage:    GC_content.pl -i fasta -o output(ID\tlength\tGC%)
Function: Calculate GC content of fasta, and ouput each slid bin GC content
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -b bin size, default 10 kb
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-22
Notes:    
\n/
    )
}

