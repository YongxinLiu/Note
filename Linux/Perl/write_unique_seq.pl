#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'l:o:d:h:', \%opts );
&usage unless ( exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{l}=6 unless defined($opts{l});
print "Input length is $opts{l}\nOutput file is $opts{o}\n";

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
open OUTPUT,">$opts{o}";
my @base=('A','C','G','T');
foreach  (@base) {
	$a=$_;
	foreach	(@base){
		$b=$_;
		foreach	(@base){
			$c=$_;
			foreach	(@base){
				$d=$_;
				foreach	(@base){
					$e=$_;
					foreach	(@base){
						foreach	(@base){
							print OUTPUT $a,$b,$c,$d,$e,$_,"\n";
						}
					}
				}
			}
		}
	}
}
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
Usage:    write_unique_seq.pl -l seq_length -o output_file
Function: given length to write unique barcode sequence, such AGCG
Command:  -l default is 6
          -o output file name (Must)
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-11-11
Notes:    
\n/
    )
}