#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:s:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Species is $opts{s}\n" if defined($opts{s});
$opts{s}="ath" unless defined($opts{s});
$opts{d}="Rep13" unless defined($opts{d});


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#    SW   perc perc perc  query         position in query              matching         repeat                 position in repeat
# score   div. del. ins.  sequence4      begin    end          (left)   repeat           class/family       begin   end    (left)     ID
#282     14.5    0.0     7.2     1       1       105     (30427566)      C       ATREP18 DNA     (1142)  649     561     1       
open OUTPUT,">$opts{o}";
while (<INPUT>) {
#0Chr		1		2	3start	4end	5		6chain	7		8type $1,name $2
#Chr1    TAIR10  gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Des=NAFYB
	chomp;
	$id++;
	my @tmp=split/\s+/;
	$tmp[8]="-" if $tmp[8] eq "C";
	print OUTPUT "$tmp[4]\t$opts{d}\trepeat\t$tmp[5]\t$tmp[6]\t$tmp[0]\t$tmp[8]\t.\tID=$opts{s}-rep$id;Note=$tmp[10];Annotation=$tmp[9]\n";
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
sub usage {#awk 'NR>3' all.con.out.rm2012 |sed "s/^ *//g"|sed "s/  */\t/g">osa7_rm2012.tsv
    die(
        qq/
Usage:    format_repeat_gff.pl -i inpute_file -o output_file
Function: format repeatmasker tsv result to gff, use standard ID genome RM result, or no header and table sperated format RM result
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -s species, default is ath
          -d database version, default is Rep13
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-01-13
Notes:    
\n/
    )
}