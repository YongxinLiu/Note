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
#ID0									flag1	chr2			position3	mapQ4
#HWI-ST833:261:3:2101:7972:43821#0       153     C135125790      61      50
my %database; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($database{$tmp[1]})) {
		$database{$tmp[0]}=$tmp[4] if $tmp[4]>$database{$tmp[0]};
	}else{
		$database{$tmp[0]}=$tmp[4];
	}
}
close DATABASE;


###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
open OUTPUTA,">$opts{o}A";
open OUTPUTB,">$opts{o}B";
open OUTPUTC,">$opts{o}C";
my %input;
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($input{$tmp[1]})) {
		$input{$tmp[0]}=$tmp[4] if $tmp[4]>$input{$tmp[0]};
	}else{
		$input{$tmp[0]}=$tmp[4];
	}
	
}
close INPUT;
foreach $id (keys %input) {
	if (!defined($database{$id})) { #如果B中不存在，则为A特异ID
		print OUTPUTA "$id\n";
	}elsif ($input{$id}==$database{$id}){ #如果A B质量相同，为common ID
		print OUTPUTC "$id\n";
		delete $database{$id};
	}elsif ($input{$id}>$database{$id}){ #如果A>B,则更可能来自A
		print OUTPUTA "$id\n";
		delete $database{$id};
	}elsif ($input{$id}<$database{$id}){ #如果A<B,则更可能来自B
		print OUTPUTB "$id\n";
		delete $database{$id};
	}
}

foreach $id (keys %database) { #剩余B中的ID全为B特异
	print OUTPUTB "$id\n";
}
close OUTPUTA;
close OUTPUTB;
close OUTPUTC;

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
Usage:    samID_group_mapQ.pl -i Sam1(no header) -d Sam2 -o outputA B common -h header num
Function: Read two sam file, compare each by mapQ, and sort into A B specifc and common(A B C)
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