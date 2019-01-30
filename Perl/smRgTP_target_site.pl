#!/usr/bin/perl
use strict;
use warnings;

my $miRNA;
my $annotation;
my $target;
open RESULT,">>target_site.txt";
while(<>){
	chomp;
	if (/>/) {
		my @tmp=split('\t',$_); 
		$target=(split(/\|/,$tmp[2]))[0];
		$miRNA=$tmp[0];
		$miRNA=~s/>//;
		$annotation=(split(/\|/,$tmp[2]))[2];
		$annotation=~s/unspliced-genomic //;
		next;
	}
	if (/Sbjct/) {
		m/\D+(\d+)\D+(\d+)/;
		print RESULT $target,"\t",$1,"\t",$2,"\t",$miRNA,"\t",$annotation,"\n";
	}
}

close RESULT;

#my $lastsplicingID="";
#open RESULT,">>osa_intron_site.txt";
#while (<>){
#	#LOC_Os12g43300.1|13112.m04564|intron_3  LOC_Os12g41620|13112.t03849|unspliced-genomic   91.94   248     18      2       159     406     5228    4983    4e-63    244
#	chomp;
#	my @tmp=split('\t',$_);
#	#print $tmp[0],"\n";
#	my $splicingID=(split('\|',$tmp[0]))[0];
#	#print $splicingID,"\n";
#	my $splicinglocID=(split('\.',$splicingID))[0];
#	#print $splicinglocID,"\n";
#	my $targetlocID=(split('\|',$tmp[1]))[0];
#	#print $targetlocID,"\n";
#	if ($splicinglocID eq $targetlocID) {
#		#print $targetlocID,"\n";
#		if ($lastsplicingID ne $splicingID) {
#			print RESULT "\n$splicingID\t";
#			$lastsplicingID=$splicingID;
#		}
#		print RESULT "$tmp[8]\t$tmp[9]\t";
#	}else{
#		next;
#	}
#}
#close RESULT;
