#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

##database in array
#my (@tmp1,@tmp2);
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
open OUTPUT,">$opts{o}";
#ID	length	GCnumber	GCcontent
my (%database,%total,%each);
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
print OUTPUT "ChrID\tLength\tA\tC\tG\tT\tGC\%\n";
foreach  (sort keys %database) {
#	print $database{$_},"\n";
	$gc=&gc($database{$_});
	$len=length($database{$_});
	$per=$gc/$len*100;
	printf OUTPUT "$_\t$len\t$each{A}\t$each{C}\t$each{G}\t$each{T}\t%.2f\n",$per;
}
$total=$total{A}+$total{G}+$total{C}+$total{T}+$total{N};
$gc=$total{G}+$total{C};
$per=$gc/$total*100;
printf OUTPUT "Total\t$total\t$total{A}\t$total{C}\t$total{G}\t$total{T}\t$gc\t%.2f\n",$per;
printf OUTPUT "Total\t$total\tA %.3f C %.3f G %.3f T %.3f\n",$total{A}/$total,,$total{C}/$total,,$total{G}/$total,,$total{T}/$total;
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
Function: Calculate GC content of fasta, and ouput each sequence and total GC content
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-22
Notes:    
\n/
    )
}

sub gc{
	my $seq=$_[0];
#	print $seq,"\n";
	my @nt=split//,$seq;
	%each=('A'=>0,'C'=>0,'G'=>0,'T'=>0);
	foreach  (@nt) {
#		print $_,"\n";
		$total{$_}++;
		$each{$_}++;
	}
	return $each{G}+$each{C};
}