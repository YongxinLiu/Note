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
my %list;
while (<DATABASE>) {
	next unless />/;
	chomp;
	$list{$_}=0;
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
#>cel-let-7
#UACACUGUGGAUCCGGUGAGGUAGUAGGUUGUAUAGUUUGGAAUAUUACCACCGGUGAAC
open INPUT,"<$opts{i}";
my %database;
while (<INPUT>) {
        chomp;
        if (/>/) {
			chomp;
            $chr=$_;
        }else{
            $database{$chr}.=$_;
        }
}
open OUTPUT,">$opts{o}";
foreach  $id(keys %database) {
	next if defined($list{$id});
	print OUTPUT "$id\n$database{$id}\n";
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
Usage:    fasta_subtraction.pl -i all_fasta_seq -d subtract_fasta_seq -o selected_fasta_seq
Function: select fasta seq from a database by ID, ID must same in hash
Command:  -i fasta (Must)
          -o output selected fasta (Must)
          -d list, ID in first line
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2016-11-19
Notes:    
\n/
    )
}
