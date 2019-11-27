#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Bio::FeatureIO;

my $inFile = shift;
my ($name, $path, $suffix) = fileparse($inFile, qr/\.gff/);
my $outFile = $path . $name . ".gtf";

my $inGFF = Bio::FeatureIO->new( '-file' => "$inFile",
								 '-format' => 'GFF',
								 '-version' => 3 );
my $outGTF = Bio::FeatureIO->new( '-file' => ">$outFile",
								  '-format' => 'GFF',
								  '-version' => 2.5);

while (my $feature = $inGFF->next_feature() ) {
	if ($feature->type->name eq 'exon' || $feature->type->name eq 'CDS') {
		my $parent = ($feature->get_Annotations('Parent'))[0];
		my $transcript = $parent->value;
		my ($gene) = $transcript =~ m/(\w+)\.\d+$/;
		my $transcript_id = Bio::Annotation::SimpleValue->new( '-value' => $transcript,
															   '-tagname' => 'transcript_id');
		my $gene_id = Bio::Annotation::SimpleValue->new( '-value' => $gene,
														 '-tagname' => 'gene_id');
		$feature->add_Annotation($transcript_id);
		$feature->add_Annotation($gene_id);
	}
	$outGTF->write_feature($feature);
}
