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

###############################################################################
#Read the database in memory(opt)
###############################################################################
open DATABASE,"<$opts{d}";
#Model_name      Type				Short_description
#AT1G01010.1     protein_coding  NAC domain containing protein 1
my (%tair,%cover);
while (<DATABASE>) {
#Chr0	source1	2type	start3	end4		5score	strand6	7CDS_phased	detail8(id,name,biotype)
#Chr1    omicslab        sRNA_cluster    44570   44956   26.76   -       .       ID=SRC5;length_type=24-nt
#Chr9	source10 type11	start12	end13	score14	strand15 CDS_phased16	detail17(id,name,biotype)
#Chr1    TAIR10  gene    44677   44787   .       +       .       ID=AT1G01073;Note=protein_coding_gene;Name=AT1G01073    111
	chomp;
	my @tmp=split/\t/;

	$tmp[8]=~/ID=([^;]+)/;
	$id=$1;
#	print "$tmp[18]\n";
	$tmp[17]=~/ID=([^;]+)/;
	$gi=$1;
	if (defined($tair{$id})){
		next if $cover{$id}<$tmp[18];
		$cover{$id}=$tmp[18];
		$tair{$id}=$gi;
	}else{
		$cover{$id}=$tmp[18];
		$tair{$id}=$gi;
	}
}
close DATABASE;

###############################################################################
#Main text.
###############################################################################
open INPUT,"<$opts{i}";
#est_id locus   FPKM_NL FPKM_SL FPKM_NR FPKM_SR P_NL_SL P_NR_SR P_NL_NR P_SL_SR annotation      biological_process
#XLOC_000001     Chr1:3630-5899  5.00205 5.27591 23.5308 38.8806 0.0336365696776894      6.03355050428681e-19    2.59784606653356e-123   1.29806438263668e-208   NAC domain containing protein 1;        multicellular organismal development, DNA-dependent, transport, developmental processes, transcription, regulation of transcription; 
open OUTPUT,">$opts{o}";
#TAIRID	est_id locus   FPKM_NL FPKM_SL FPKM_NR FPKM_SR P_NL_SL P_NR_SR P_NL_NR P_SL_SR annotation      biological_process
#AT1G00000	XLOC_000001     Chr1:3630-5899  5.00205 5.27591 23.5308 38.8806 0.0336365696776894      6.03355050428681e-19    2.59784606653356e-123   1.29806438263668e-208   NAC domain containing protein 1;        multicellular organismal development, DNA-dependent, transport, developmental processes, transcription, regulation of transcription; 
while ($opts{h}>0) { #filter header
	$a=<INPUT>;
	print OUTPUT "TAIRID\tCover\t$a";
	$opts{h}--;
}
my $i=0;
while (<INPUT>) {
	chomp;
	my @tmp=split/\t/;
	$tmp[1]=~/\:(\d+)\-(\d+)/;
	$len=$2-$1+1;
	if (defined($tair{$tmp[0]})) {
		$coverage=$cover{$tmp[0]}/$len*100;
		printf OUTPUT "$tair{$tmp[0]}\t%.0f\t$_\n",$coverage;
	}else{
		$i++;
		printf OUTPUT "ATXG%05d\t0\t$_\n",$i;
	}
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
Usage:    intersect_gff_RNAseqtair.pl -i gff_cluster  -d intersect -o result
Function: Show the intersect file as a more comprehensive read file,Reads maize gene annotation from phytozome zma annotation, change transcription ID to gene ID; reads each cluster annotation in array, only keep the biggest overlap; result each cluster in one line.
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
		  -h header
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2013-08-23
Notes:    
\n/
    )
}