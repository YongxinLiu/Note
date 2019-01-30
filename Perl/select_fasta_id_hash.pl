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
#>cel-let-7
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
foreach  $list(@list) {
	next unless defined($database{$list});
	print OUTPUT ">$list\n$database{$list}\n";
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
Usage:    select_fasta_id_hash.pl -i all_fasta_seq -d list -o selected_fasta_seq
Function: select fasta seq from a database by ID, ID must same in hash
Command:  -i fasta (Must)
          -o output selected fasta (Must)
          -d list, ID in first line
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2015-03-24
Notes:    
\n/
    )
}
