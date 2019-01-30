#!/usr/bin/perl
use strict;

if ($ARGV[0] =~ /^-[h?]/) {
print "Usage: determine-phred FILE

Reads a sam or fastq, possibly gzipped and returns the phred-scale,
either 64 or 33, based on a quick scan of the data in the file.
";
exit 0;
}
my $cnt;
my $dphred = 64;
if ($ARGV[0] =~ /\.gz$/) {
$ARGV[0] = "gunzip -c '$ARGV[0]'|";
}
my $qual;
my $comm;
my $fmt;
if (@ARGV > 1) {
my @mult = @ARGV;
for my $f (@mult) {
@ARGV = ($f);
determine();
print "$f\t$dphred\n";
}
} else {
determine();
print "$dphred\n";
}

sub determine {
$_ = <>;
if (/^\@/ && ! /^\@SQ\t/) {
# fastq
scalar <>; # read
$comm = scalar <>; # comment
if (!(substr($comm,0,1) eq '+')) {
die "Unknown file format\n";
}
$qual = <>;
chomp $qual;
$fmt = 'fq';
} else {
# sam
$fmt = 'sam';
$qual = (split(/\t/, $_))[10];
}
if (!$qual) {
die "Unknown file format\n";
}

my $ssiz=7000; # sample size
while($qual) {
for (my $i =length($qual)/2; $i < length($qual); ++$i) {
if (ord(substr($qual,$i,1)) < 64) {
$dphred = 33;
$cnt=$ssiz; # last
last;
}
}
$qual = '';
last if ++$cnt >= $ssiz; # got enough
if ($fmt eq 'fq') {
# fastq
scalar <>; # id
scalar <>; # read
scalar <>; # comment
$qual = <>;
chomp $qual;
} else {
# sam
$qual = (split(/\t/, $_))[10];
}
}
}
