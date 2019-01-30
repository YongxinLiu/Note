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
#open DATABASE,"<$opts{d}";
#my %database; #database in hash
#while (<DATABASE>) {
#	chomp;
#	my @tmp=split/\t/;
#	$database{$tmp[1]}=$tmp[2];
#}

#my (@tmp1,@tmp2); #database in array
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
#>S000494589 uncultured bacterium; YRM60L1D06060904	Lineage=Root;rootrank;Bacteria;domain;"Actinobacteria";phylum;Actinobacteria;class;Acidimicrobidae;subclass;Acidimicrobiales;order;"Acidimicrobineae";suborder;Acidimicrobiaceae;family;Acidimicrobium;genus$
open OUTPUT,">$opts{o}_seq.fa";
#>S000494589
open OUTPUT2,">$opts{o}_tax.txt";
# S000494589 k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__
while ($opts{h}>0) { #filter header
	<INPUT>;
	$opts{h}--;
}
while (<INPUT>) {
#	chomp;
	if (/>/) {
		my %tax;
		my @tmp=split/\t/; # tab split before Lineage
		@tmp1=split(/\s+/,$tmp[0]); # split first ID and last species
		print OUTPUT $tmp1[0],"\n";
#		print $tmp1[-1],"\n";
		$tmp1[0]=~s/>//g;
		$tmp[1]=~/(\w+)\;domain/;
		$tax{domain}=$1;
		$tmp[1]=~/(\w+)\"\;phylum/;
		$tax{phylum}=$1;
		$tmp[1]=~/(\w+)\;class/;
		$tax{class}=$1;
		$tmp[1]=~/(\w+)\;order/;
		$tax{order}=$1;
		$tmp[1]=~/(\w+)\;family/;
		$tax{family}=$1;
		$tmp[1]=~/(\w+)\;genus/;
		$tax{genus}=$1;
		print OUTPUT2 "$tmp1[0]\tk__$tax{domain}; p__$tax{phylum}; c__$tax{class}; o__$tax{order}; f__$tax{family}; g__$tax{genus}; s__$tmp1[-1]\n";
		

	}else{
		print OUTPUT uc($_);
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
Usage:    format_rdp2qiime.pl -i current_Bacteria_unaligned.fa -o current_Bacteria_unaligned.txt.fasta -h header num
Function: format rdp current_Bacteria_unaligned fasta file into qiime used fasta & taxonomy
Command:  -i inpute file name (Must)
          -o output file name (Must)
          -d database file name
          -h header line number, default 0
Author:   Liu Yong-Xin, woodcorpse\@163.com, QQ:42789409
Version:  v1.0
Update:   2017-03-21
Notes:    
\n/
    )
}