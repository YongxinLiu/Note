#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:k:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}"; # fasta file
open OUTPUT,">$opts{o}";
#chr0	start1	end2	copies3
#2       20319180        20319200        3
open INPUT,"<$opts{i}";
my @chr; # chr list
my %db; # genome database
my %seq; # short reads
my %copy; # copy number
open OUTPUT,">$opts{o}";
while (<INPUT>) {
	chomp;
	if (/>/) {	# read chr list into @chr; seq into %db
		$_=~/>([^\s]+)/;
		$chr=$1;
		push @chr,$1;
	}else{
#       $_=~s/U/T/g;
		$db{$chr}.=$_;
	}
}
foreach  $chr(@chr) { # statistic copy number
	$start=0;
	$len=length($db{$chr});
	while ($start+$opts{k}<=$len) {
		$seq{$chr}[$start]=substr($db{$chr},$start,$opts{k});
		$seq=$seq{$chr}[$start];
		$copy{$seq}++;
		$seq=~tr/ACGT/TGCA/;
		$seq=reverse($seq);
		$copy{$seq}++;
		$start++;
	}
}

foreach  $chr(@chr) { # output position, copies, which include GG
	$start=0;
	$len=length($db{$chr});
	while ($start+$opts{k}<=$len) {
		$seq{$chr}[$start]=substr($db{$chr},$start,$opts{k});
		$seq=$seq{$chr}[$start];
		next if $copy{$seq}>3;
		next unless ($seq=~/^CC/ || $seq=~/GG$/);
		print OUTPUT "$chr\t$start\t",$start+$opts{k},"\t$copy{$seq}\n";
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
Usage:    sequence_copies.pl -i genome.fa -o postion_count.bedgraph -k 23
Function: calculate each K bp sequence copies
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -k reads length
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-03-24
Notes:    
\n/
    )
}