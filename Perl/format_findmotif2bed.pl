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
$opts{h}=1 unless defined($opts{h});

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
#0FASTA ID        1Offset  2Sequence        3Motif Name      4Strand  5MotifScore
#9_56173200_56173400     -50     CAAAAAAACAGAAA  8-TTTCTGTTTTTTTV        +       15.639300
#9_56173100_56173300     50      CAAAAAAACAGAAA  8-TTTCTGTTTTTTTV        +       15.639300
#9_56173000_56173200     -62     ATTGTGTTTTTTTC  8-TTTCTGTTTTTTTV        -       14.071328
#9_56172900_56173100     38      ATTGTGTTTTTTTC  8-TTTCTGTTTTTTTV        -       14.071328
open OUTPUT,">$opts{o}";
#0chr	1start	2end	3name		4reads	5strand	
#9       989616  989862  CAAAAAAACAGAAA    15.639300       +   
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	my @pos=split/\_/,$tmp[0];
	if ($tmp[4] eq '+') {
		$start=($pos[1]+$pos[2])/2+$tmp[1];
		$end=$start+length($tmp[2]);
	}else{
		$end=($pos[1]+$pos[2])/2+$tmp[1]+1;
		$start=$end-length($tmp[2])+1;
	}
	print OUTPUT "$pos[0]\t$start\t$end\t$tmp[2]\t$tmp[5]\t$tmp[4]\n";
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
Usage:    format_findmotif2bed.pl -i find_motif_result -o bed_file -h 1
Function: format findmotif result, and change position to bed
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 1, remove title line
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    find_motif_result ID must only seq position such as "9:56173200-56173400" (create by write_bed_bin.pl and fastaFromBed), Offset is distance to middle of seq, bed file start position is less 1 than real.
\n/
    )
}