#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:t:', \%opts );
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
my %cultured; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$cultured{$tmp[0]}="YES";
}
close DATABASE;

open TAX,"<$opts{t}";
my %tax; #database in hash
while (<TAX>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[1]=~s/\.//g; # 删除所有.，graphlan有点会错误分隔
	$tmp[1]=~s/\w__//g; # 删除所有级别
	$tmp[1]=~s/;/\t/g; # 按制表符分割
	$tax{$tmp[0]}=$tmp[1];
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
print OUTPUT "OTU_ID\trecovered\tkindom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n";
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	print OUTPUT "$tmp[0]\t";
	if (defined($cultured{$tmp[0]})) {
		print OUTPUT "YES\t";
	}else{
		print OUTPUT "NO\t";
	}
	print OUTPUT $tax{$tmp[0]},"\n";
	
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
        qq!
Usage:    graphlan_culture.pl -i OTU_ID -o output_file -d cultured_list -t taxonomy -h header num
Function: Add cultured and taxonomy in select ID fro graphlan
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017/10/30
Notes:    
\n!
    )
}