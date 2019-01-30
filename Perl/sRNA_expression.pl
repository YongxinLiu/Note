#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:c:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{c}=1 unless defined($opts{c});


###############################################################################
#Main text.
###############################################################################
my $i=0;
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	$header[$i]=<INPUT>;
	chomp($header[$i]);
	$opts{h}--;
}
my $hash_sRNA; #record all sRNA and expression
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
#srid    seq     loci    len
#Sr0007491       AGATGACGTTGTTGACTTTTTCTC        5       24
	$hash_sRNA{$tmp[1]}[$i]=$_;
}

close INPUT;
my @filelist=glob "$opts{d}";
foreach $file(@filelist){
	my %list;
	my $total=0;
	$i++;
	$j=$i+1;
	open DATABASE,"<$file";
	$header[$i]=basename($file);
	while (<DATABASE>) { #record one sample in hash
		chomp;
		my @tmp=split/\t/;
		$list{$tmp[0]}=$tmp[1];
		$total+=$tmp[1];
	}
	$header[$j]=$total;
	foreach  $sRNA(keys %hash_sRNA) {
		if (defined($list{$sRNA})) {
			$hash_sRNA{$sRNA}[$i]=$list{$sRNA};
			$hash_sRNA{$sRNA}[$j]=$list{$sRNA}*1000000/$total;
		}else{
			$hash_sRNA{$sRNA}[$i]=0;
			$hash_sRNA{$sRNA}[$j]=0;
		}
	}
	$i=$j;
	close DATABASE;
}
for (0..$i) {
	print OUTPUT $header[$_],"\t";
}
print OUTPUT "\n";
	
foreach  $sRNA(keys %hash_sRNA) {
	for (0..$i) {
		print OUTPUT $hash_sRNA{$sRNA}[$_],"\t";
	}
	print OUTPUT "\n";
}

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
Usage:    sRNA_expression.pl -i inpute_file -o output_file -d database	-h header num
Function: calculate sRNA reads and RPM, based on sRNA list and sRNA database. Rule, read sRNA list in hash-array, read each sample and then add reads and RPM to hash-arry. Using a header save head info.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -c column select, default 1
          -h header number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2014-05-13
Notes:    
\n/
    )
}