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
open DATABASE,"<$opts{d}";#reads fasta in hash
#>AT1G51370.2 | Symbols:  | F-box/RNI-like/FBD-like domains-containing protein | chr1:19045615-19046748 FORWARD LENGTH=346
#MVGGKKKTKICDKVSHEEDRISQLPEPLISEILFHLSTKDSVRTSALSTKWRYLWQSVPGLDLDPYASSNTNTIVSFVES
#FFDSHRDSWIRKLRLDLGYHHDKYDLMSWIDAATTRRIQHLDVHCFHDNKIPLSIYTCTTLVHLRLRWAVLTNPEFVSLP
my %database;
my $chr;
while (<DATABASE>) {
	chomp;
	if (/>/) {
		$chr=$_;
#		print $chr,"\n";
	}else{
		$database{$chr}.=$_;
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";

foreach  $id(keys %database) {
	if ($database{$id}=~/R\w{2}L\w{4}N/) {
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
Usage:    fetch_protein.pl -i inpute_file -o output_file -d database-h header num
Function: give a motif or seq, and find all the protein in database, X represent each AA
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2013-08-22
Notes:    
\n/
    )
}