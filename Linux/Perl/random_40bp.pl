#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;


###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));


###############################################################################
#Read the genome in hash
###############################################################################
open DATABASE,"<$opts{d}";
my @database;
my @chromosome;
my $n=-1;
while (<DATABASE>) {
	chomp;
	if (/>/) {
		$n++;
		$database[$n]="";
		$chromosome[$n]=(split('\s',$_))[0];
	}else{
		$database[$n].=$_;
	}
}
close DATABASE;

my @seqlength;
my $total=0;
my @percent;
foreach (@database){
	push @seqlength,length($_);
	$total+=length($_);
}
my $l=0;
foreach (@database){
	$l+=length($_);
	push @percent,$l/$total;
}


###############################################################################
#Main text.
###############################################################################
open OUTPUT,">$opts{o}";
for (1..$opts{i}){
	my $start=rand();
	for (0..@percent){
		if ($start<$percent[$_]){
			fetchseq($_,$percent[$_]-$start);
			last;
		}
	}
}

close OUTPUT;


sub fetchseq{
	my ($id,$percent)=@_;
	my $start=int($percent*$seqlength[$id]-40);
	if($start<0){
		$start=0;
	}
	my $seq=substr($database[$id],$start,40);
	print OUTPUT "$chromosome[$id]\t$start\n$seq\n";
}

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
Usage:    template.pl -i inpute_file -o output_file
Function:
Command:  -i inpute file name (Must)
          -o output file name
          -d database file name
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2011-09-20
Notes:    
\n/
    )
}