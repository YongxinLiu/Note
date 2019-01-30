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
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[1];
#	push @tmp2,@tmp[2];
#}
#close DATABASE;

#open a list file
#my %list;
#my @filelist=glob "$opts{i}";
#foreach $file(@filelist){
#	open DATABASE,"<$file";
#	$file=basename($file);
#	while (<DATABASE>) {
#		my @tmp=split/\t/;
#		$list{$file}{nr}++;
#	}
#	close DATABASE;
#}

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
my %sgRNA;
while (<INPUT>) {
#0Qname	1Flag	2Rname	3Pos	4MapQ	5Cigar	6Rnext	7Pnext	8Tlen	9Seq	10Qual	11mismatch	12mismatch_site	13mismatch&gap
#zb7-sg2 16      1       273004611       255     19M     *       0       0       GCTGATACATAGCATCTTG     IIIIIIIIIIIIIIIIIII     XA:i:0  MD:Z:19 NM:i:0
	chomp;
	next if /^@/;
	my @tmp=split/\t/;
	next unless defined($tmp[13]);
	@tmp1=split(/:/,$tmp[13]);
	#print $tmp1[2],"\n";
	$sgRNA{$tmp[0]}[$tmp1[2]]++;
}

print OUTPUT "ID\tOff-target0\tOff-target1\tOff-target2\n";
foreach  (sort keys %sgRNA) {
	$sgRNA{$_}[0]-=1;
	$sgRNA{$_}[1]=0 unless defined($sgRNA{$_}[1]);
	$sgRNA{$_}[2]=0 unless defined($sgRNA{$_}[2]);
	print OUTPUT "$_\t$sgRNA{$_}[0]\t$sgRNA{$_}[1]\t$sgRNA{$_}[2]\n";
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
Usage:    sgRNA_off-target.pl -i inpute_file.sam -o output_file.txt
Function: show each gene off-target sites
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