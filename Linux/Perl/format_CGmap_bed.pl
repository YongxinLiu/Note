#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:t:c:p:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{t}="C" unless defined($opts{t});
$opts{p}=0 unless defined($opts{p});



###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;

#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr0	strand1	position2		type3			methylation5
#3       G       218921469       CHH     CC      0.333333333333  2       6
#3       G       218921471       CG      CG      1.0     6       6
#3       G       218921472       CHG     CC      0.333333333333  2       6
open OUTPUT,">$opts{o}";
#chr0		start1			end2		id3		methylation4	strand5
#3       218922760       218922760       CHG27   0.428571429     -
#3       218923200       218923200       CHG28   0.6     +
#3       218923372       218923372       CHG29   1       -
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my $i=0;
my $strand;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[3]=~/$opts{t}/; # select 5mC type 
	$chr=$tmp[0];
	$chr=$opts{c} if defined($opts{c}); # rename chr ID
	$start=$tmp[2]-$opts{p}; # change position
	$i++;
	if ($tmp[1] eq "C") {
		$strand="+";
	}else{
		$strand="-";
	}
	print OUTPUT "$opts{c}\t$start\t$start\t$tmp[3]$i\t$tmp[5]\t$strand\n";
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
Usage:    format_CGmap_bed.pl -i inpute_file -o output_file -d database	-h header num
Function: fomrat bs_seeker CGmap to bed
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
          -t type CG, CHG, CHH, default "C"
          -c chrID, default not change
          -p posoition change, default -
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-10-28
Notes:    
\n/
    )
}