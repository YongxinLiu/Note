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
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my %seq;
while (<INPUT>) {
	chomp;
	next if />/;
	@tmp=split(//,$_);
	for $i(0..$#tmp) {
#		print "$i\t$tmp[$i]\n";
		$seq{$i}{$tmp[$i]}++;
	}
}

foreach $i(0..$#tmp) {
	$seq{$i}{A}=0 unless defined($seq{$i}{A});
	$seq{$i}{G}=0 unless defined($seq{$i}{G});
	$seq{$i}{C}=0 unless defined($seq{$i}{C});
	$seq{$i}{T}=0 unless defined($seq{$i}{T});
	if ($seq{$i}{A}>=$seq{$i}{G} && $seq{$i}{A}>=$seq{$i}{C} && $seq{$i}{A}>=$seq{$i}{T}) {
		printf OUTPUT "A\t%.2f\n",$seq{$i}{A}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		print "A";
		#print "A\t$seq{$i}{A}\t",($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100,"\t",$seq{$i}{A}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100,"\n";
		next;
	}elsif ($seq{$i}{G}>=$seq{$i}{A} && $seq{$i}{G}>=$seq{$i}{C} && $seq{$i}{G}>=$seq{$i}{T}) {
		printf OUTPUT "G\t%.2f\n",$seq{$i}{G}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		print "G";
		#printf "G\t%.2f\n",$seq{$i}{G}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		next;
	}elsif ($seq{$i}{C}>=$seq{$i}{A} && $seq{$i}{C}>=$seq{$i}{G} && $seq{$i}{C}>=$seq{$i}{T}) {
		printf OUTPUT "C\t%.2f\n",$seq{$i}{C}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		print "C";
		#printf "C\t%.2f\n",$seq{$i}{C}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		next;
	}elsif ($seq{$i}{T}>=$seq{$i}{A} && $seq{$i}{T}>=$seq{$i}{C} && $seq{$i}{T}>=$seq{$i}{G}) {
		printf OUTPUT "T\t%.2f\n",$seq{$i}{T}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		print "T";
		#printf "T\t%.2f\n",$seq{$i}{T}/($seq{$i}{A}+$seq{$i}{G}+$seq{$i}{C}+$seq{$i}{T})*100;
		next;
	}
}
print "\n";


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
Usage:    template.pl -i inpute_file -o output_file -d database	-h header num
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-10
Notes:    
\n/
    )
}