#!/usr/bin/perl -w
use POSIX qw(strftime);
use Getopt::Std;
use File::Basename;

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my %opts;
getopts( 'i:o:d:h:r:g:c:', \%opts );
&usage unless ( exists $opts{i} && exists $opts{o} );
my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is $opts{i}\nOutput file is $opts{o}\n";
print "Database file is $opts{d}\n" if defined($opts{d});
$opts{r}=1000000 unless defined($opts{r});
$opts{g}="Zma3" unless defined($opts{g});

###############################################################################
#Read the database in memory(opt)
###############################################################################
#open DATABASE,"<$opts{d}";
#my (@tmp1,@tmp2);
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	push @tmp1,$tmp[0];
#}
#close DATABASE;



###############################################################################
#Main text.
###############################################################################
#select bam chromosome, and format to bed
#my $bed_list;
#foreach $chr(@tmp1) {
#	`bamtools convert -in $opts{i} -region $chr -format bed -out $chr\.bed`;
#	$bed_list.="$chr\.bed ";
#}
#merge each chr, and delete temp
#`cat $bed_list>$opts{o}\.bed`;
#`rm $bed_list`;
#$scale=1000000/$opts{r};
#`bedtools genomecov -i $opts{o}\.bed -g $opts{d} -scale $scale -bg>$opts{o}\.wig`;
#`rm $opts{o}\.bed`;

#bamtools filter -in sDic15.bam -region 9 -out chr9_sDic15.bam # mainly study chr5, filter chr5
#bedtools genomecov -ibam chr9_sDic15.bam -scale 0.0378744983049268 -bg>bedgraph/chr9_sDic15.bedgraph
#~/software/IGVTools/igvtools toTDF bedgraph/chr9_sDic15.bedGraph tdf/chr9_sDic15 Zma3


#$dir=dirname($opts{o});
#print "mkdir $dir\n";
#`mkdir $dir`;
#print "bamtools filter -in $opts{i} -region $opts{c} -out $opts{o}\.bam\n";
#`bamtools filter -in $opts{i} -region $opts{c} -out $opts{o}\.bam`;
$scale=1000000/$opts{r};
print "bedtools genomecov -i $opts{i} -scale $scale -g $opts{h} -bg>$opts{o}\.bedgraph\n";
`bedtools genomecov -i $opts{i} -scale $scale -g $opts{h} -bg>$opts{o}\.bedgraph`;
print "~/software/IGVTools/igvtools toTDF $opts{o}\.bedgraph $opts{o} $opts{g}\n";
`~/software/IGVTools/igvtools toTDF $opts{o}\.bedgraph $opts{o} $opts{g}`;
print "rm $opts{o}\.bam $opts{o}\.bedgraph\n";
`rm $opts{o}\.bedgraph`;


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
Usage:    bam_select_chr_tdf.pl -i inpute_file_bam -o output_file_tdf -r mapped_total_reads -c chromsome_ID -g genome_name
Function: format bam to tdf, and select chr parameter c
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -r sample size or mapped reads, use for scale to RPM, default=1000000, no scale
          -c chromosome number
		  -g genome name, such tair10, Zma3
		  -h genome each chr length
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2014-08-01
Notes:    
\n/
    )
}