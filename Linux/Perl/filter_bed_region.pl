#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:p:s:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{s}=1 unless defined($opts{s});

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{p}";
#chr	start	end	length	(unit in bp)
#1	134124000	135251000	301476924
my %db; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$db{$tmp[0]}{start}=$tmp[1]/$opts{s};
	$db{$tmp[0]}{end}=$tmp[2]/$opts{s};
}

close DATABASE;


###############################################################################
#Main text.
###############################################################################
my(@hang,$site,%num,$Bd,$n,$start,$end,$m);
%num=();
open F1,"<$opts{i}";   ### bed-format
#1	184549146	184549546	Nucleosome:694475	0.0258537
#1	38797081	38797481	Nucleosome:145474	0.0430896
#1	53477181	53477581	Nucleosome:200822	0.0344717
open OUT,">$opts{o}";
 while(<F1>){
	my @tmp=split/\t/;
	$m=$tmp[0];
	$n=$tmp[1];
	print "$m\t$n\n";
	next if $n<$db{$m}{start} || $n>$db{$m}{end};
	print OUT "$_";
 }
close OUT;

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
Usage:    filter_bed_region.pl -i gff -o gff -p centromere_position
Function: filter centromere region annotation, also can filter any regions annotation
Command:  -i inpute file name (Must)
          -o output file name (Must)
		  -p position include chr start end chr_length
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-06-11
Notes:    
\n/
    )
}