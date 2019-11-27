#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:s:t:m:e:', \%opts );
&usage unless ( exists $opts{i});
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\n";
$opts{h}=0 unless defined($opts{h});
$opts{s}="./" unless defined($opts{s});
$opts{t}="./" unless defined($opts{t});
$opts{m}="mv" unless defined($opts{m});
$opts{e}="fna" unless defined($opts{e});


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#SampleID	OldID
#CviPollen	BN2-886
#ColStigma1	BN2-900
while ($opts{h}>0) { #filter header
        <INPUT>;
        $opts{h}--;
}
while (<INPUT>) {
	chomp;
	@tmp=split/\t/;
	if ($opts{m} eq "mv") {
	`$opts{m} $opts{s}$tmp[1]_1.fq.gz $opts{t}$tmp[0]_1.fq.gz`;
	print "$opts{m} $opts{s}$tmp[1]_1.fq.gz $opts{t}$tmp[0]_1.fq.gz\n";
	}elsif ($opts{m} eq "ln") {
#	`$opts{m} $opts{s}$tmp[1]/$tmp[1].extendedFrags.fastq $opts{t}$tmp[0].fq`;
#	print "$opts{m} $opts{s}$tmp[1]/$tmp[1].extendedFrags.fastq $opts{t}$tmp[0].fq \n";
	`$opts{m} $opts{s}$tmp[1]/$tmp[1].$opts{e} $opts{t}$tmp[0].$opts{e}`;
	`sed 's/$tmp[1]/$tmp[0]/' $opts{t}$tmp[0].$opts{e}`;
	print "$opts{m} $opts{s}$tmp[1]/$tmp[1].$opts{e}  $opts{t}$tmp[0].$opts{e}  \n";
	}elsif ($opts{m} eq "cp") {
	`$opts{m} $opts{s}$tmp[1]/$tmp[1].$opts{e} $opts{t}$tmp[0].$opts{e}`;
	`sed -i "s/$tmp[1]/$tmp[0]/" $opts{t}$tmp[0].$opts{e}`;
	print "$opts{m} $opts{s}$tmp[1]/$tmp[1].$opts{e}  $opts{t}$tmp[0].$opts{e}  \n";
	print "sed \"s/$tmp[1]/$tmp[0]/\" $opts{t}$tmp[0].$opts{e}\n";
	}else{
		print "Nothing done!!!";
	}
}
close INPUT;

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
        qq!
Usage:    rename_fq.pl -i inpute_list
Function: mv or link all column2(OldID) to column1(NewID) according input list for single fastq files
Command:  -i inpute file name (Must)
          -e extended, default fna
          -s source directory, default current
          -t target directory, default current
          -m method, default mv, alternative ln
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2018/3/17
Notes:    
\n!
    )
}