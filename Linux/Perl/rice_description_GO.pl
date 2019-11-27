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
open DATABASE,"<$opts{d}";
my %db; #database in hash
#loci0	GO1	type2	GO_name3	4	5
#LOC_Os01g01010.1        GO:0030234      F       enzyme regulator activity       IEA     TAIR:AT3G59570
#LOC_Os01g01010.1        GO:0007165      P       signal transduction     IEA     TAIR:AT3G59570
#LOC_Os01g01010.1        GO:0006139      P       nucleobase, nucleoside, nucleotide and nucleic acid metabolic process   IEA     TAIR:AT3G59570
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$db{$tmp[0]}{ID}.="$tmp[1]\;";
	$db{$tmp[0]}{type}.="$tmp[2]\;";
	$db{$tmp[0]}{name}.="$tmp[3]\;";
}


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#chr0	locus1	model2	start3	stop4	ori5	is_TE6	is_expressed7	is_representative8	annotation9
#Chr1	LOC_Os01g01010	LOC_Os01g01010.1	2903	10817	+	N	Y	Y	TBC domain containing protein, expressed
#Chr1	LOC_Os01g01010	LOC_Os01g01010.2	2984	10562	+	N	Y	N	TBC domain containing protein, expressed
#Chr1	LOC_Os01g01019	LOC_Os01g01019.1	11218	12435	+	N	Y	Y	expressed protein
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$a=$tmp[2];
	$db{$a}{ID}="" unless defined($db{$a}{ID});
	$db{$a}{type}="" unless defined($db{$a}{type});
	$db{$a}{name}="" unless defined($db{$a}{name});
	print OUTPUT "$a\t$tmp[9]\t$db{$a}{ID}\t$db{$a}{type}\t$db{$a}{name}\n";
	
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
Usage:    rice_description_GO.pl -i all.locus_brief_info.7.0 -o GI_descript_GOID_GOtype_GOname -d all.GOSlim_assignment	-h 1 header num
Function: add rice GO term 
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