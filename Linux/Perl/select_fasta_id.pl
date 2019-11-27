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
###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
my @list;
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	push @list,$tmp[0];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
#>cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
#UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
open INPUT,"<$opts{i}";
my %database;
open OUTPUT,">$opts{o}";
while (<INPUT>) {
        chomp;
        if (/>/) {
            $_=~/>([^\s]+)/;
            $chr=$1;
#           print $chr,"\n";
        }else{
			$_=~s/U/T/g;
            $database{$chr}.=$_;
        }
}
foreach  $id(keys %database) {
	foreach  $list(@list) {
		if ($id eq $list) {
			print OUTPUT ">$id\n$database{$id}\n";
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
Usage:    select_fasta_id.pl -i all_fasta_seq -d list -o selected_fasta_seq
Function: select fasta seq from a database by ID
Command:  -i fasta (Must)
          -o output selected fasta (Must)
          -d list, new ID in one line
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-07-20
Notes:    
\n/
    )
}
