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
# ID0	ID1	Evalue2	Score3
# TRIUR3_26600-T1 EMT16955        0.0      1269
my %database; #database in hash
while (<DATABASE>) {
	chomp;
	my @tmp=split/\t/;
	$database{$tmp[0]}{besthit}=$tmp[1]; #��ID��Ӧ��ID�ʹ��
	$database{$tmp[0]}{score}=$tmp[3];
	if (defined($database{$tmp[1]})) { #��ID�Ի�������ߴ�ֵĶ�
		if ($tmp[3]>$database{$tmp[1]}{score}) {
			$database{$tmp[1]}{besthit}=$tmp[0];
			$database{$tmp[1]}{score}=$tmp[3];
		}
	}else{
		$database{$tmp[1]}{besthit}=$tmp[0];
		$database{$tmp[1]}{score}=$tmp[3];
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
# ID0	ID1	Evalue2	Score3
# TRIUR3_26600-T1 EMT16955        0.0      1269
open OUTPUT,">$opts{o}";
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if (defined($database{$tmp[0]})) {
		if ($tmp[3]>$database{$tmp[0]}{score}) {#ֻ�޸ĸ����ʵĶ�
			$database{$tmp[0]}{besthit}=$tmp[1];
			$database{$tmp[0]}{score}=$tmp[3];
			if (defined($database{$tmp[1]})) {
				if ($tmp[3]>$database{$tmp[1]}{score}) {
					$database{$tmp[1]}{besthit}=$tmp[0];
					$database{$tmp[1]}{score}=$tmp[3];
				}
			}else{
				$database{$tmp[1]}{besthit}=$tmp[0];
				$database{$tmp[1]}{score}=$tmp[3];
			}
		}
	}else{
		$database{$tmp[0]}{besthit}=$tmp[1];
		$database{$tmp[0]}{score}=$tmp[3];
		if (defined($database{$tmp[1]})) {
			if ($tmp[3]>$database{$tmp[1]}{score}) {
				$database{$tmp[1]}{besthit}=$tmp[0];
				$database{$tmp[1]}{score}=$tmp[3];
			}
		}else{
			$database{$tmp[1]}{besthit}=$tmp[0];
			$database{$tmp[1]}{score}=$tmp[3];
		}
	}
}
close INPUT;
foreach  $id(keys %database) {
	next if !defined($database{$id});
	my $idB=$database{$id}{besthit};
	next if !defined($database{$idB});
	if ($database{$idB}{besthit} eq $id) {
		print OUTPUT "$id\t$idB\t$database{$id}{score}\n$idB\t$id\t$database{$idB}{score}\n";
		delete $database{$id};
		delete $database{$idB};
	}
}
open OUTPUTT,">$opts{o}tmp";
foreach  $id(keys %database) {
	next if !defined($database{$id});
	my $idB=$database{$id}{besthit};
	next if !defined($database{$idB});
	print OUTPUTT "$id\t$idB\t$database{$id}{score}\n$idB\t$database{$idB}{besthit}\t$database{$idB}{score}\n";
	delete $database{$id};
	delete $database{$idB};
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
Usage:    best_hit_mutual.pl -i first file -o output file -d second file -h header num
Function: get the mutual best hit gene pair, use for find homolog
		  output is 1:1 homolog pair, output_tmp include none 1:1 homolog
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-09-11
Notes:    
\n/
    )
}