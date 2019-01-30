#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:c:l:m:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
$opts{c}=0 unless defined($opts{c});
$opts{l}=36 unless defined($opts{l});
$opts{m}=16 unless defined($opts{m});

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
my ($total_seq, $total_reads, $select_seq, $select_reads) =(0,0,0,0);
while (<INPUT>){
	chomp;
	$total_seq++;
	my ($seq,$count)=split/\t/;
	next unless $count=~/\d+/;
	$total_reads+=$count;
	next unless $count>$opts{c};
	next unless length($seq)<=$opts{l};
	next unless length($seq)>=$opts{m};
	print OUTPUT "$total_seq&$count\t$seq\n";
	$select_seq++;
	$select_reads+=$count;
}

print "Total sequences $total_seq\nSelect sequence $select_seq\nTotal reads $total_reads\nSelect reads $select_reads\n";

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
Usage:    format_sRNA_psmap.pl -i inpute_table_reads_count -o output_file
Function: format non-redundancy sRNA from fasta to psmap (format psRobot_map)
Command:  -i tab file (Must)
          -o psRobot map file
          -c count filter, default 0, not filter
          -l length max filter, default 36, recommend 26
          -m length min filter, default 16, recommend 18
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2014-06-02
Notes:    Add read and length filter
\n/
    )
}
