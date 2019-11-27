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
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
##database in hash
#my %database;
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

##database in array
#my (@tmp1,@tmp2);
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
#match0	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#     	match	match	   	count	bases	count	bases	      	name9     	size10	start	end	name13     	size14	start	end	count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#31	0	0	0	0	0	0	0	+	ath-SRC1	201	0	31	ath-SRC22801	101	20	51	1	31,	0,	20,
open OUTPUT,">$opts{o}";
my (%len,%overlap);
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$len{$tmp[9]}=$tmp[10];
	$len{$tmp[13]}=$tmp[14];
#	if (defined($overlap{$tmp[9]}{$tmp[13]})) {
#		if ($overlap{$tmp[9]}{$tmp[13]}>$tmp[0]) {
#			next;
#		}
#	}
	$overlap{$tmp[9]}{$tmp[13]}=$tmp[0];
	$overlap{$tmp[13]}{$tmp[9]}=$tmp[0];
	
}
my @src;
@src=sort keys %overlap;
print "Totally ", $#src+1," SRCs\n";
foreach  (@src) {
#	print $keya,"\n";
#	next if $hash{$keya}{len}==0;
#	print OUTPUT "$keya\t$hash{$keya}{len}";
#	$hash{$keya}{len}=0;
#	foreach  $keyb(sort keys %hash) {
#		next unless defined($hash{$keyb}{len});
#		next if $hash{$keyb}{len}==0;
#		if (defined($hash{$keya}{$keyb})) {
#			print OUTPUT "\t$keyb\t$hash{$keyb}{len}\t$hash{$keya}{$keyb}";
#			$hash{$keyb}{len}=0;
#		}
#		next if $hash{$keyb}{len}==0;
#		if (defined($hash{$keyb}{$keya})) {
#			print OUTPUT "\t$keyb\t$hash{$keyb}{len}\t$hash{$keyb}{$keya}";
#			$hash{$keyb}{len}=0;
#		}
#	}
	next if $len{$_}==0;
	&cluster($_);
	print OUTPUT "\n";

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
Usage:    cluster_blat.pl -i inpute_file -o output_file
Function: Template for Perl
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-05-22
Notes:    
\n/
    )
}


sub cluster{
	$key1=$_[0];
	my @key2=sort keys %{$overlap{$key1}};
#		print "$key1\t";
#		foreach (@key2) {
#			print "$_\t";
#		}
#		print "\n";
	foreach $key2(@key2) {
#		print "$key1\t$key2\n";
		&printo($key1,$key2);
		delete $overlap{$key2}{$key1};
		delete $overlap{$key1}{$key2};
		&cluster($key2);
	}
}

sub printo{
	@tmp=@_;
	if ($len{$tmp[0]}>0) {
		print OUTPUT "$tmp[0]\t";
		$len{$tmp[0]}=0;
	}
	if ($len{$tmp[1]}>0) {
		print OUTPUT "$tmp[1]\t";
		$len{$tmp[1]}=0;
	}
}