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
open OUTPUT,">$opts{o}";

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

my $miRNA="";
my $start="";
my $gene="";
while (<INPUT>) {
	chomp;
	my @input=split('\t',$_);
	if (($miRNA eq $input[0]) && $start==$input[1]) {
		next;
	}else{
		$miRNA=$input[0];
		$start=$input[1];
		($gene)=split('\.',$input[4]);
#		print $gene,"\n";
	}
#	print $_,"\n";
	my $count=0;
	my @loc;
	foreach my $database(@database) {
#		print $database,"\n";
		if ($database=~/$gene/) {
			$count++;
			push @loc,$database;
		}
	}
	if ($count<2) {
		next;
	}
	my $overlap=0;
	foreach my $loc(@loc) {
		my @intron_site=split('\t',$loc);
		shift @intron_site;
		foreach my $site (@intron_site) {
			if ($input[1]> $site && $input[2]<$site) {
				$overlap++;
			}
		}
	}
	if ($overlap<$count) {
		print OUTPUT $_,"\n";
		foreach my $loc(@loc) {
			my $unoverlap=2;
			my @intron_site=split('\t',$loc);
			shift @intron_site;
			foreach my $site (@intron_site) {
				if ($input[1]> $site && $input[2]<$site) {
					$unoverlap=1;
				}
			}
			my $demension=scalar @intron_site;
#			print $demension,"\n";
#			print @intron_site,"\n";
			for (my $i=0;$i<$demension ;$i=$i+2) {
#				print $intron_site[$i],"\t",$input[2],"\t",$intron_site[$i+1],"\t",$input[1],"\n";
				if ($intron_site[$i]<$input[2] && $intron_site[$i+1]>$input[1]) {
					$unoverlap=3;
				}
			}
			#isoform1 quantify¡úor¡ý, proportion ¡ú¡ý.
			#isoform2 quantify¡ý, proportion¡ý.
			#isoform3 quantify¡ú, proportion¡ú¡ü.

			if ($unoverlap==1) {
				print OUTPUT "=-\t",$loc,"\n";
			}elsif ($unoverlap==3){
				print OUTPUT "=+\t",$loc,"\n";
			}else{
				print OUTPUT "--\t",$loc,"\n";
			}
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
Usage:    template.pl -i inpute_file -o output_file
Function: Filter miRNA potential regulate alternative splicing region
		  Add the alternative splicing classify.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2011-10-10
Notes:    
\n/
    )
}