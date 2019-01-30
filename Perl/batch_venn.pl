#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:w:u:E:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{d} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{h}=0 unless defined($opts{h});
$opts{i}='result_k1-c/group_venn.txt' unless defined($opts{i});
$opts{d}='result_k1-c/family.txt' unless defined($opts{d});
$opts{w}=4 unless defined($opts{w});
$opts{u}=4 unless defined($opts{u});

###############################################################################
#Main text.
###############################################################################
if (-e "$opts{i}") {
open INPUT,"<$opts{i}";
# filter header
 while ($opts{h}>0) {
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	if ($#tmp==1) {
		print "sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1]\n";
		`sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1]`;
		print "sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1]\n";
		`sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1]`;
		print "vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1]\n";
		`vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1]`;
	}elsif ($#tmp==2) {
		print "sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2]\n";
		`sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2]`;
		print "sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2]\n";
		`sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2]`;
		print "vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1],$tmp[2]\n";
		`vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1],$tmp[2]`;
	}elsif ($#tmp==3) {
		print "sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3]\n";
		`sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3]`;
		print "sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3]\n";
		`sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3]`;
		print "vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1],$tmp[2],$tmp[3]\n";
		`vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1],$tmp[2],$tmp[3]`;
	}elsif ($#tmp>3) {
		print "sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3] -g $tmp[4]\n";
		`sp_vennDiagram.sh -E pdf -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3] -g $tmp[4]`;
		print "sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3] -g $tmp[4]\n";
		`sp_vennDiagram.sh -E png -w $opts{w} -u $opts{u} -f $opts{d} -a $tmp[0] -b $tmp[1] -c $tmp[2] -d $tmp[3] -g $tmp[4]`;
		print "vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4]\n";
		`vennNumGenerator.py -i $opts{d} -g $tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4]`;
	}	
}
close INPUT;
}
###############################################################################
#Record the program running time!
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";

###############################################################################
#Scripts usage
###############################################################################
sub usage {
    die(
        qq!
Usage:    batch_venn.pl -i venn_overlap_group_list -o work_dir -d DA_gene_or_otu_list -h header num
Function: batch call sp_vennDiagram.sh, according to group_venn.txt, also report each overlap number by vennNumGenerator.py
Command:  -i venn_overlap_group_list (Must)
          -o work_dir, same with -d work dir
          -d DA_gene_or_otu_list  (Must)
          -h header line number, default 0
          -w width of figure, default 2
          -u height of figure, default 2
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.1
Update:   2017/6/5
Notes:    1.0 draw venn diagram
		1.1 Add vennNumGenerator.py report overlap list, add venn legnth and height
\n!
    )
}