#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:c:s:e:p:n:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{p}=0 unless defined($opts{p});
###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#3       Rep13   repeat  1       2943    21583   -       .       ID=ath-rep1;Note=LTR/Copia;Annotation=Copia13-ZM_I-int
#3       Rep13   repeat  2944    4283    10056   -       .       ID=ath-rep2;Note=LTR/Copia;Annotation=Copia13-ZM_LTR
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	next unless $tmp[0] eq $opts{c}; # select chr
	next unless $tmp[3]>$opts{s}; # discard < start
	last if $tmp[3]>$opts{e}; # stop in end
	next if $tmp[4]>$opts{e};
	if (defined($opts{n})) {
		print OUTPUT "$opts{n}\t$tmp[1]\t$tmp[2]\t",$tmp[3]-$opts{p},"\t",$tmp[4]-$opts{p},"\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[8]\n" 
	}else{
		print OUTPUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t",$tmp[3]-$opts{p},"\t",$tmp[4]-$opts{p},"\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[8]\n" 
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
Usage:    select_gff.pl -i gff -o part_gff -c chromosome -s start -e end -h header_num -n chr_new_name -p reposition_coordinate
Function: select part file from big gff, also can rename chr and change coordinate
Command:  -i gff file (Must)
          -o gff file (Must)
		  -c chr ID (Must)
		  -s start (Must)
		  -e end (Must)
          -n new chr name, default no change
		  -p default=0, no change
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-10-15
Notes:    
\n/
    )
}