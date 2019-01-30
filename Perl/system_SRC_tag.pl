#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:', \%opts );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));

##############################################################################
#Main text.
###############################################################################
`SRC_tag.pl -i bio/SRCID.00 -d difference/sig_Pol4_flower.diff0 -o bio/SRCID.01 -t PolIV-flower`;
`SRC_tag.pl -i bio/SRCID.01 -d difference/sig_Pol4_seedling.diff0 -o bio/SRCID.02 -t PolIV-seedling`;
`SRC_tag.pl -i bio/SRCID.02 -d difference/sig_Pol5_flower.diff0 -o bio/SRCID.03 -t PolV-flower`;
`SRC_tag.pl -i bio/SRCID.03 -d difference/sig_Pol5_seedling.diff0 -o bio/SRCID.04 -t PolV-seedling`;
`SRC_tag.pl -i bio/SRCID.04 -d difference/sig_RDR1_leaf.diff0 -o bio/SRCID.05 -t RDR1-leaf`;
`SRC_tag.pl -i bio/SRCID.05 -d difference/sig_RDR2_flower.diff0 -o bio/SRCID.06 -t RDR2-flower`;
`SRC_tag.pl -i bio/SRCID.06 -d difference/sig_RDR2_seedling.diff0 -o bio/SRCID.07 -t RDR2-seedling`;
`SRC_tag.pl -i bio/SRCID.07 -d difference/sig_RDR6_flower.diff0 -o bio/SRCID.08 -t RDR6-flower`;
`SRC_tag.pl -i bio/SRCID.08 -d difference/sig_RDR6_leaf.diff0 -o bio/SRCID.09 -t RDR6-leaf`;
`SRC_tag.pl -i bio/SRCID.09 -d difference/sig_RDR6_seedling.diff0 -o bio/SRCID.10 -t RDR6-seedling`;
`SRC_tag.pl -i bio/SRCID.10 -d difference/sig_DCL1_leaf.diff0 -o bio/SRCID.11 -t DCL1-leaf`;
`SRC_tag.pl -i bio/SRCID.11 -d difference/sig_DCL2_seedling.diff0 -o bio/SRCID.12 -t DCL2-seedling`;
`SRC_tag.pl -i bio/SRCID.12 -d difference/sig_DCL3_seedling.diff0 -o bio/SRCID.13 -t DCL3-seedling`;
`SRC_tag.pl -i bio/SRCID.13 -d difference/sig_DCL4_seedling.diff0 -o bio/SRCID.14 -t DCL4-seedling`;
`SRC_tag.pl -i bio/SRCID.14 -d difference/sig_IP_AGO1_flower.diff0 -o bio/SRCID.15 -t AGO1-flower`;
`SRC_tag.pl -i bio/SRCID.15 -d difference/sig_IP_AGO1_leaf.diff0 -o bio/SRCID.16 -t AGO1-leaf`;
`SRC_tag.pl -i bio/SRCID.16 -d difference/sig_IP_AGO1_root.diff0 -o bio/SRCID.17 -t AGO1-root`;
`SRC_tag.pl -i bio/SRCID.17 -d difference/sig_IP_AGO2_seedling.diff0 -o bio/SRCID.18 -t AGO2-seedling`;
`SRC_tag.pl -i bio/SRCID.18 -d difference/sig_IP_AGO4_flower.diff0 -o bio/SRCID.19 -t AGO4-flower`;
`SRC_tag.pl -i bio/SRCID.19 -d difference/sig_IP_AGO4_leaf.diff0 -o bio/SRCID.20 -t AGO4-leaf`;
`SRC_tag.pl -i bio/SRCID.20 -d difference/sig_IP_AGO4_root.diff0 -o bio/SRCID.21 -t AGO4-root`;
`SRC_tag.pl -i bio/SRCID.21 -d difference/sig_IP_AGO4_seedling.diff0 -o bio/SRCID.22 -t AGO4-seedling`;
`SRC_tag.pl -i bio/SRCID.22 -d difference/sig_IP_AGO9_flower.diff0 -o bio/SRCID.23 -t AGO9-flower`;
`SRC_tag.pl -i bio/SRCID.23 -d difference/sig_IP_AGO10_flower.diff0 -o bio/SRCID.24 -t AGO10-flower`;
`cp bio/SRCID.24 bio/SRC_FC4_RPM5_P1Q5.id #SRC annotation, fold change > 2, P-value <0.01, FDR <0.05`;

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
Usage:    system.pl -i inpute_file -o output_file
Function: batch run command
Command:  -i inpute file name (Must)
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-06-03
Notes:    
\n/
    )
}