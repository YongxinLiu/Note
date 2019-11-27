#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Bio::SeqIO;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
open INPUT,"<$opts{i}";
open OUTPUT,">>$opts{o}";

###############################################################################
#Main text.
###############################################################################
open DATABASE,"<$opts{d}";
my @database;
while (<DATABASE>) {
	chomp;
	push @database,$_;
#	print $_,"\n";
}
close DATABASE;

while (<INPUT>) {
	chomp;
	my @input=split('\t',$_);
#	print $_,"\n";
	foreach my $database(@database) {
#		print $database,"\n";
		if ($database=~/$input[0]/) {
#			print $input[0],"\n";
			my @intron_site=split('\t',$database);
			shift @intron_site;
#			my $alternative_splicing=shift @intron_site;
#			print $alternative_splicing,"\n";
			foreach my $site (@intron_site) {
#				print $site,"\t",$input[1],"\t",$input[2],"\n";
				if ($input[1]> $site && $input[2]<$site) {
#					$input[4]=~s/unspliced-genomic //;
					print OUTPUT "$input[3]\t$input[1]\t$input[2]\t$input[4]\t$database\n";
#					last;
				}
			}
		}
	}
}

###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
close INPUT;
close OUTPUT;

###############################################################################
#Scripts usage and about.
###############################################################################
sub usage {
    die(
        qq/
Usage:    filter_overlap_branchpoint.pl -i inpute_file -o output_file -d database_file
Function: Filter the miRNA target region overlap with mRNA transcript branch point.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name (Must)
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-09-20
Notes:    Input format:LOC_Os01g02470	493	472	osa-miR169g(f)	unspliced-genomic conserved hypothetical protein
          Database format:LOC_Os01g01010.1	367	451	715	1454	1554	2554	2659	4233	5043	5125	5249	5329	5419	5505	5707	6307	6716	7201	7286	7371	7529	7601	
          Output format:
\n/
    )
}