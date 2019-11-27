#!/usr/bin/perl
# 29 September 2014
# If using this pipeline please cite : XXXXXXXXXX
#-----------------------------------------------------------+
#                                                           |
# TASR - Transposon Annotation using Small RNAs.            |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: MOAINE EL BAIDOURI                               |
# CONTACT: moaine.elbaidouri@gmail.com                      |        
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
# VERSION: V.1.0                                            |
#                                                           |                                               
#-----------------------------------------------------------+

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/shuffle/;
use Statistics::Descriptive;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::PrimarySeq;
use Bio::SeqIO;

## setup defaults options ##

my $ref;
my $siRNAs_file;
my $usearchv;
my $blastv       = 'NCBI-BLAST+';
my $out_folder   = 'TASR_output';
my $temp	 = 'delete';
my $siRNAs       = 4;
my $CPU          = 1;
my $minlen       = 80;
my $maxlen       = 20000;
my $idenblast	 = 70;
my $evalue	 = 1e-20;
my $window       = 150;
my $idensilix    = 0.8;
my $overlapsilix = 0.8;
my $copy_num     = 2;
my $help         = 0;

GetOptions(
    'ref=s'          => \$ref,
    'sfile=s'        => \$siRNAs_file,
    'outfold=s'      => \$out_folder,
    'blastv=s'       => \$blastv,
    'usearchv=s'     => \$usearchv,
    'nsirna=i'       => \$siRNAs,
    'win=i'          => \$window,
    'minlen=i'       => \$minlen,
    'maxlen=i'       => \$maxlen,
    'idenblast=i'    => \$idenblast,
    'evalue=s'       => \$evalue,
    'cpu=i'          => \$CPU,
    'idclust=s'      => \$idensilix,
    'overlapclust=s' => \$overlapsilix,
    'cnumber=i'      => \$copy_num,
    'temp=s'	     => \$temp,
    'help!'          => \$help,
) or die "Incorrect usage!\n";

if ($help) {
    print "

TASR version v.1.0

usage: perl TASR_v.1.0.pl -ref genome.fasta -sfile siRNAs.fasta -usearchv path/usearch7 [options]



	-ref          (mandatory)
	 		   	  reference genome sequence in fasta format
		
	-sfile        (mandatory)
			  	  siRNAs sequences in fasta format

        -usearchv     (mandatory)
			 	  please specify the path to executable file name (ex : -usearchv /usr/bin/usearch7.0.1090_i86linux32)
	
	-outfold      (optional)
			 	  output directory (default : -outfold TASR_output)

	-cpu          (optional)
				  minimum number of threads used for mapping and blast (default : -cpu 1)

	-nsirna       (optional)
			 	  minimum number of siRNAs to consider TEs candidates loci (default : -nsirna 4)

	-win          (optional)
				  maximum window size in bp to merge mapped siRNAs (default : -win 150)

	-minlen       (optional)
				  minimum interval length in bp to consider TEs candidates loci default : -minlen 80)
	
	-maxlen       (optional)
				  maximum interval length in bp to consider TEs candidates loci (default : -maxlen 20000)

	-idclust      (optional)
				  minimum % identity to accept blast hits for building families (in [0,1]) (default : -idclust 0.8)

	-overlapclust (optional)
				  minimum % overlap to accept blast hits for building families (in [0,1]) (default : -overlapclust 0.8)

	-cnumber      (optional)
				  minimum copy number to consider TEs families (default : -cnumber 2 )

	-idenblast    (optional)
				  minimum identity percentage cut-off for blast (default : -idenblast 70)

	-evalue       (optional)
				  Expected value (default : -evalue 1e-20)

	-Blastv       (optional)
			 	  If using older blast version please specify '-Blastv NCBI-BLAST' (default : -Blastv NCBI-BLAST+)

	-temp	      (optional)
	 		          temporary files : keep or delete (default : -temp delete)
      
	-h --help
				  print this menu
";
}
if ( !$ref || !$siRNAs_file || !$usearchv ) {
    print
        "\n\nWARNING : reference genome, siRNAs sequences and usearch file name are required. Please type : 'perl TASR_v.1.0.pl -h' for help \n";
    exit 1;
}
##########################################################################
##################### Get siRNAs mapping intervals #######################
##########################################################################

### Test program and perl modules ####

my $BL=qx(which blast2);
my $BLplus=qx(which blastn);
if ($BL eq "" && $BLplus eq ""){
    print "\n\nPlease install blast2 or blast+ before runing TASR\\n\n\n";
    exit;
}
my $BWT=qx(which bowtie2);
if ($BWT eq ""){
    print "\n\nPlease install bowtie2 before runing TASR\n\n\n";
    exit;
}
my $TRF=qx(which trf);
if ($TRF eq ""){
    print "\n\nPlease install Tandem Repeat Finder before runing TASR\n\n\n";
    exit;
}
my $SILIX=qx(which silix);
if ($SILIX eq ""){
    print "\n\nPlease install silix before runing TASR\n\n\n";
    exit;
}
my $USEARCH=qx(ls $usearchv);
if ($USEARCH eq ""){
    print "\n\nNo usearch executable file found ! please verify the correct path to uesearch program\n\n\n";
    exit;
}

## log folder ##
if ( -d "TASR_log" ) { system("rm -r TASR_log"); }
system("mkdir TASR_log");

## Print options ###

open( OPTIONS, ">TASR_log/TASR_options" );
print OPTIONS " -reference genome     (-ref) 		: $ref		\n
  -siRNAs file          (-sfile)     	: $siRNAs_file	\n 	     
    
  -number of sirnas     (-nsirna) 	: $siRNAs   	\n
  -windows size	        (-win)		: $window       \n 
  -min length           (-minlen)	: $minlen       \n
  -max length           (-maxlen)	: $maxlen  	\n
  -clustering identity  (-idclust)      : $idensilix	\n
  -min overlap          (-overlapclust) : $overlapsilix	\n
  -min copy number      (-cnumber)	: $copy_num	\n
  -min Blast identity   (-idenblast)	: $idenblast	\n
  -Blast evalue 	(-evalue)	: $evalue	\n

  -output folder        (-outfold)	: $out_folder	\n 
  -number of threads    (-cpu)		: $CPU		\n
  -usearch version      (-usearchv)	: $usearchv   	\n
  -Blast version        (-Blastv)    	: $blastv   	\n\n\n\n";

## reference genome ##

RENAME($ref);
my $db = Bio::DB::Fasta->new("$ref.rename");
my %header;

## Create temporary directory ##

my $time = time;
system("mkdir tmp$time");

## Mapping siRNAs against the reference genome ###
print "\n** Bowtie2 genome index ** \n";
system("bowtie2-build $ref.rename reference >/dev/null 2>/dev/null");

print "\n        ..done.. \n";

system("mv reference*.bt2 tmp$time");

print "\n   ** siRNAs mapping **\n";

system(
    "bowtie2 -k 400 -D 10 -R 5 -N 1 -L 15 -i S,1,0.50 -p $CPU -x tmp$time/reference -f $siRNAs_file -S tmp$time/24.sirna.sam  2>TASR_log/bowtie2.log"
);

print "\n        ..done..  \n";

## Sam to bam to bed and sorting ##

SAMTOBED("tmp$time/24.sirna.sam");
system(
    "sort -k1,1 -k2,2n tmp$time/24.sirna.sam.bed >tmp$time/24.sirna.sam.sorted.bed" 
); #--parallel=$CPU 

## Merging siRNAs intervals using window ##

MERGEBED( "tmp$time/24.sirna.sam.sorted.bed", "$window" );

## Filter intervals using minlen,maxlen and nsiran ####

open( MERGE, ">tmp$time/24.sirna.merged.$window.$siRNAs.filtered.bed" );
open( LEN,   "<tmp$time/24.sirna.sam.sorted.bed.merged" );
while ( my $lLEN = <LEN>) {
    chomp($lLEN);
    my @lLEN  = split( " ", $lLEN );
    my $sinum = $lLEN[3];
    my $inlen = $lLEN[2] - $lLEN[1];
    if (   ( $sinum >= $siRNAs )
        && ( $inlen >= $minlen )
        && ( $inlen <= $maxlen ) )
    {
        print MERGE "$lLEN\n";
    }
}
close(LEN);
close(MERGE);

### Extract fasta from intervals ###

BEDTOFASTA( "tmp$time/24.sirna.merged.$window.$siRNAs.filtered.bed",
    "$ref.rename" );

## Mask simple repeat using TRF program

print "\n** Masking simple repeats **\n";

system(
    "trf tmp$time/24.sirna.merged.$window.$siRNAs.filtered.bed.fasta 2 7 7 80 10 50 500 -m -h >/dev/null 2>/dev/null"
);
system(
    "rm 24.sirna.merged.$window.$siRNAs.filtered.bed.fasta.2.7.7.80.10.50.500.dat"
);
system(
    "mv 24.sirna.merged.$window.$siRNAs.filtered.bed.fasta.2.7.7.80.10.50.500.mask tmp$time/24.sirna.merged.$window.$siRNAs.mask.fasta"
);

print "\n         ..done..  \n";

# refilter sequence length after masking #

system("sed -i 's.N..g' tmp$time/24.sirna.merged.$window.$siRNAs.mask.fasta");
open( FILTER, ">tmp$time/24.sirna.merged.$window.$siRNAs.mask.final" );

my $seq_in = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => "tmp$time/24.sirna.merged.$window.$siRNAs.mask.fasta"
);
while ( my $seq1 = $seq_in->next_seq() ) {
    my $id1 = $seq1->primary_id;
    chomp $id1;
    my $seqrelen = $seq1->seq;
    chomp $seqrelen;
    my $lseq = length($seqrelen);
    if ( $lseq >= $minlen && $lseq <= $maxlen ) {
        print FILTER ">", $id1, "\n", $seqrelen, "\n";
    }
}
close(FILTER);

## reformat trf output for silix ##

my $seqout = Bio::SeqIO->new(
    -file     => ">tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta",
    '-format' => 'Fasta'
);
my $seqin = Bio::SeqIO->new(
    -file     => "tmp$time/24.sirna.merged.$window.$siRNAs.mask.final",
    '-format' => 'Fasta'
);
while ( my $seqref = $seqin->next_seq ) {
    $seqout->write_seq($seqref);
}

## Blast all against all ##

print "\n** Blast all against all **\n";


if ( $blastv eq 'NCBI-BLAST' ) {
    system(
        "formatdb -i tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta -p F >TASR_log/formatdb.blast.log"
    );
    system(
        "blast2 -a $CPU -p blastn -d tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta -i tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta -o tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta.blast -m 8 -r 2 -e $evalue -P $idenblast -F F  >TASR_log/blast.log"
    );
}
else {
    system(
        "makeblastdb -in tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta -dbtype nucl >TASR_log/formatdb.blast.log"
    );
    system(
        "blastn -num_threads $CPU -db tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta -query tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta -out tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta.blast -outfmt 6 -perc_identity $idenblast -evalue $evalue  >TASR_log/blast.log"
    );
}
print "\n        ..done..  \n";

## clustering ##

print "\n** Clustering families **\n";

system(
    "silix -i $idensilix -r $overlapsilix tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta.blast -f FAM > tmp$time/clusters 2>TASR_log/silix.log"
);
system("mkdir tmp$time/clustering");
system(
    "silix-split tmp$time/24.sirna.merged.$window.$siRNAs.mask.final.fasta tmp$time/clusters -o tmp$time/clustering/ -n $copy_num >/dev/null 2>/dev/null"
);
print "\n 	..done..\n";

##########################################################################
##################### Redefine TEs boundaries #############################
##########################################################################

print "\n ** Redefining TE boundaries **\n";

my @id_last;
my $id_last;
my @pos;
my $hit_name;
my $query_name;
my @files = glob("tmp$time/clustering/*.fasta");
foreach (@files) {
    chomp($_);
    my @name = split( "\/", $_ );
    my $name = pop(@name);

##### Calculate median length of the TEs #####

    my @len = qx(fastalength tmp$time/clustering/$name);
    my @leng;
    foreach (@len) {
        chomp($_);
        my @ll = split( " ", $_ );
        push( @leng, $ll[0] );
    }

    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@leng);
    my $quartile4 = $stat->quantile(3);

############ Extend left and right $med+medplus and write new fasta file ##############

    my $outextend = Bio::SeqIO->new(
        -file   => ">tmp$time/clustering/$name.tfa",
        -format => "fasta"
    );
    my $seq_in = Bio::SeqIO->new( -format => 'fasta', -file => $_ );
    while ( my $seq = $seq_in->next_seq() ) {
        my $id = $seq->primary_id;
        chomp $id;
        my $sequence = $seq->seq;
        chomp $seq;
        my @lll = split( "_", $id );
        my $chr = $lll[0];
        if ( $quartile4 > 2000 ) {
            $quartile4 = 2000;
        }
        my $star = int( $lll[1] - ( 5 * $quartile4 ) );
        if ( $star < 0 ) {
            $star = 1;
        }
        my $en = int( $lll[2] + ( 5 * $quartile4 ) );
        my $sequen = $db->subseq( $chr, $star, $en );
        my $newid = join( "_", $chr, $star, $en, $id );
        my $fast = Bio::PrimarySeq->new( -id => $newid, -seq => $sequen );
        $outextend->write_seq($fast);
    }

##########Detect Repeated sequences using a blast all against all######

    my $copy_number = qx(grep -c ">" $_);

    if ( $copy_number <= 10 ) {
        REPEAT("tmp$time/clustering/$name.tfa");
    }
    if ( $copy_number > 10 && $copy_number < 30 ) {
        RANDOM_FASTA( "tmp$time/clustering/$name.tfa", 10 );
        REPEAT("tmp$time/clustering/$name.tfa.subset");
        system(
            "mv tmp$time/clustering/$name.tfa.subset.final.repeat tmp$time/clustering/$name.tfa.final.repeat"
        );
    }

    if ( $copy_number >= 30 ) {
        RANDOM_FASTA( "tmp$time/clustering/$name.tfa", 20 );
        REPEAT("tmp$time/clustering/$name.tfa.subset");
        system(
            "mv tmp$time/clustering/$name.tfa.subset.final.repeat tmp$time/clustering/$name.tfa.final.repeat"
        );
    }
    system(
        "$usearchv -id 0.8 -cluster_fast  tmp$time/clustering/$name.tfa.final.repeat -centroids tmp$time/clustering/$name.tfa.final.repeat.centroid>/dev/null 2>/dev/null"
    );

    if ( $blastv eq 'NCBI-BLAST' ) {
        system(
            "formatdb -i tmp$time/clustering/$name.tfa.final.repeat.centroid -p F >TASR_log/formatdb.blast.log"
        );
        system(
            "blast2 -a $CPU -p blastn -d tmp$time/clustering/$name.tfa.final.repeat.centroid -i tmp$time/clustering/$name.tfa -o tmp$time/clustering/$name.bl -e $evalue -r 2 -F F -P $idenblast >TASR_log/blast.log"
        );
    }
    else {
        system(
            "makeblastdb -in tmp$time/clustering/$name.tfa.final.repeat.centroid  -dbtype nucl >TASR_log/formatdb.blast.log"
        );
        system(
            "blastn -num_threads $CPU -db tmp$time/clustering/$name.tfa.final.repeat.centroid -query tmp$time/clustering/$name.tfa -out tmp$time/clustering/$name.bl -perc_identity 70 -evalue $evalue >TASR_log/blast.log"
        );
    }

################################################################################

    DETECT_BORDER( "tmp$time/clustering/$name", $copy_number );





}    #foreach
print "\n	..done..\n";

if ( !$out_folder ) {
    system("rm -r TASR_output");
    system("mkdir TASR_output");
    system("mv tmp$time/clustering/*.fa TASR_output/");
}
else {
    system("mkdir $out_folder");
    system("mv tmp$time/clustering/*.fa $out_folder/");
}
if ($temp eq 'delete') {
   system("rm -r tmp$time");
}
system("rm $ref.rename*");

print "\n-->TASR has finished runing<--\n";

###############################
########## SUBROUTINES ########
###############################

########### SAM TO BED ########

sub SAMTOBED {
    my ($samfile) = @_;
    open( SAMBED, ">$samfile.bed" );
    open( SAM,    "<$samfile" );
    while (<SAM>) {
        if ( $_ =~ 'AS:i' ) {
            chomp($_);
            my @cols   = split( "\t", $_ );
            my $chr    = $cols[2];
            my $start  = $cols[3];
            my $seqlen = length( $cols[9] );
            my $end    = $start + $seqlen;
            print SAMBED "$chr\t$start\t$end\n";
        }
    }
    close(SAM);
    close(SAMBED);
}

##### MERGE BED ##############

sub MERGEBED {
    my ( my $bed, $window ) = @_;
    my $intv_ref;
    my $intv_start;
    my $intv_stop;
    my $intv_count = 0;
    open( BED,   "<$bed" );
    open( MERGE, ">$bed.merged" );
    while ( my $line = <BED> ) {
        chomp($line);

        my ( $bed_ref, $bed_start, $bed_stop ) = split( /\t/, $line );

        if ( !defined($intv_ref) || $bed_ref ne $intv_ref ) {
            if ( $intv_count > 0 ) {
                print( MERGE
                        "$intv_ref\t$intv_start\t$intv_stop\t$intv_count\n" );
            }

            $intv_ref   = $bed_ref;
            $intv_start = $bed_start;
            $intv_count = 0;
        }

        if ( $intv_count == 0 ) {
            $intv_start = $bed_start;
            $intv_count = 1;
        }

        else {
            if ( $bed_start - $intv_stop <= $window ) {
                $intv_count++;
            }

            else {
                print( MERGE
                        "$intv_ref\t$intv_start\t$intv_stop\t$intv_count\n" );

                $intv_start = $bed_start;
                $intv_count = 1;
            }
        }

        if ( $intv_count == 1 || $bed_stop > $intv_stop ) {
            $intv_stop = $bed_stop;
        }
    }
    close(BED);

    if ( $intv_count > 0 ) {
        print( MERGE "$intv_ref\t$intv_start\t$intv_stop\t$intv_count\n" );
    }
    close(MERGE);
}

##### RANDOM SELECT FASTA ####

sub RANDOM_FASTA {
    my ( $ffile, $sample ) = @_;
    my $seq_in = Bio::SeqIO->new( -format => 'fasta', -file => "$ffile" );
    my @array;
    while ( my $seqran = $seq_in->next_seq() ) {
        my $idrand = $seqran->primary_id;
        push( @array, $idrand );
    }
    my @sample = ( shuffle(@array) )[ 0 .. $sample - 1 ];
    my $dbran  = Bio::DB::Fasta->new("$ffile");
    my $outsample
        = Bio::SeqIO->new( -file => ">$ffile.subset", -format => "fasta" );
    foreach my $idsample (@sample) {
        my $seqsample = $dbran->seq("$idsample");
        my $seqoutsample
            = Bio::PrimarySeq->new( -id => $idsample, -seq => $seqsample );
        $outsample->write_seq($seqoutsample);
    }
}

##### DETECT REPEAT USING BLAST ####

sub REPEAT {
    my ($filename) = @_;
    if ( $blastv eq 'NCBI-BLAST' ) {
        system("formatdb -i $filename -p F >TASR_log/formatdb.blast.log");
        system(
            "blast2 -a $CPU -p blastn -d  $filename -i $filename -o $filename.bl -P $idenblast -m 8 -r 2 -e $evalue -F F >TASR_log/blast.log"
        );
    }
    else {
        system(
            "makeblastdb -in $filename -dbtype nucl >TASR_log/formatdb.blast.log"
        );
        system(
            "blastn -num_threads $CPU -db $filename -query $filename -out $filename.bl -outfmt 6 -perc_identity $idenblast -evalue $evalue >TASR_log/blast.log"
        );
    }

    open( OUTPUT, ">$filename.bed" );
    open( INPUT,  "<$filename.bl" );
    while ( my $line = <INPUT> ) {
        chomp($line);
        my @array_r = split( "\t", $line );
        if (   ( $array_r[0] ne $array_r[1] )
            && ( $array_r[3] >= 50 ) )
        {
            my @STARTEND = ( $array_r[6], $array_r[7] );
            my @srt = sort { $a <=> $b } (@STARTEND);
            print OUTPUT "$array_r[0]\t$srt[0]\t$srt[1]\n";
        }
    }
    close(INPUT);
    close(OUTPUT);
    system(
        "sort -k1,1 -k2,2n $filename.bed > $filename.bed.sorted"
    );# --parallel=$CPU
    BEDTOFASTA( "$filename.bed.sorted", "$filename" );
    system("mv $filename.bed.sorted.fasta $filename.final.repeat");
}

##### REDEFINE BOUNDARIES ####

sub DETECT_BORDER {
    my ( $name, $copy_number ) = @_;
    my $thresh;

##########################count hit number by repeat###################

    my @hits;
    my $in0 = new Bio::SearchIO( -format => 'blast', -file => "$name.bl" );
    while ( my $result0 = $in0->next_result ) {
        my $numhit0 = $result0->num_hits;
        if ($numhit0) {
            while ( my $hit0 = $result0->next_hit ) {
                my $hsp0 = $hit0->next_hsp;
                if ( defined $hsp0 ) {
                    my $hit_length   = $hit0->length;
                    my $hit_length30 = 0.3 * $hit_length;
                    if (   ( $hsp0->length('total') >= 50 )
                        && ( $hsp0->length('total') > $hit_length30 ) )
                    {
                        my $hit_name0 = $hit0->name;
                        push( @hits, $hit_name0 );
                    }

                }
            }
        }
    }

##### calculate the number of repeat threshold to keep for blast position###

    my @repeatnum;
    my %hash;
    my %counts;
    my $threshold;
    for (@hits) {
        $counts{$_}++;
    }
    foreach my $keys ( keys %counts ) {
        $hash{$keys} = $counts{$keys};
        push( @repeatnum, $counts{$keys} );
    }
    my $stat1 = Statistics::Descriptive::Full->new();
    $stat1->add_data(@repeatnum);
    my $quartile4_1 = $stat1->quantile(4);
    $threshold = $quartile4_1 * 0.30;

    if ( $copy_number > 10 ) {
        $thresh = $threshold;
    }
    if ( $copy_number <= 10 && $copy_number > 6 ) {
        $thresh = 4;
    }
    if ( $copy_number <= 6 && $copy_number >= 3 ) {
        $thresh = 3;
    }
    if ( $copy_number < 3 ) {
        $thresh = 2;
    }

## Open the blast file##

    my $out = Bio::SeqIO->new( -file => ">$name.fa", -format => "fasta" );
    my $DB = Bio::DB::Fasta->new("$name.tfa");
    my $blast_io
        = Bio::SearchIO->new( -file => "$name.bl", -format => 'blast' );

    while ( my $result = $blast_io->next_result ) {
        my $query_name   = $result->query_name;
        my @query        = split( "_", $query_name );
        my $query_length = $result->query_length;
        my @pos;
        open( OUTFILE, ">$name.hits" );
        while ( my $hit = $result->next_hit ) {
            $hit_name = $hit->name;
            if ( $hash{$hit_name} >= $thresh ) {
                while ( my $hsp = $hit->next_hsp ) {
                    if ( $hsp->length('total') >= 50 ) {
                        if ( $hit_name ne $query_name ) {
                            my $sb = $hsp->start('query');
                            my $eb = $hsp->end('query');
                            my @l  = ( $sb, $eb );
                            my @ls = sort { $a <=> $b } (@l);
                            my $s  = shift(@ls);
                            my $e  = pop(@ls);
                            print OUTFILE "$query_name\t$s\t$e\n";
                            push( @pos, $s, $e );
                        }
                    }
                }
            }
        }
        close(OUTFILE);
        system(
            "sort -k1,1 -k2,2n $name.hits > $name.hits.sorted"
        );# --parallel=$CPU
        MERGEBED( "$name.hits.sorted", "1" );
        my $FINAL_ID;
        my $seq;
        my $fasta;
        my $pos_number = scalar(@pos);

        if ( $pos_number < 2 ) {
            my $CHROM = $query[3];
            my $start = $query[4];
            my $end   = $query[5];
            $FINAL_ID = "$CHROM" . "_" . "$start" . "_" . "$end";
            $seq = $db->subseq( $CHROM, $start, $end );
            $fasta = Bio::PrimarySeq->new( -id => $FINAL_ID, -seq => $seq );
            $out->write_seq($fasta);
        }
        else {
            open( HITS, "$name.hits.sorted.merged" );
            my @arrayseq;
            my @positionss;
            my $i = 1;
            while ( <HITS> ) {
                my @HITS = split( "\t", $_ );
                my @positions = ( $HITS[1], $HITS[2] );
                my @orderbrut = sort { $a <=> $b } (@positions);
                my @order     = sort { $a <=> $b } (@positions);
                my $start     = shift(@order);
                my $end       = pop(@order);
                my @GPOS          = split( "_", $query_name );
                my $genomic_start = $GPOS[1] + $start - 1;
                my $genomic_end   = $GPOS[1] + $end - 1;
                my $genomic_chr   = $GPOS[0];
                $FINAL_ID
                    = "$genomic_chr" . "_"
                    . "$genomic_start" . "_"
                    . "$genomic_end";
                $seq = $DB->subseq( $query_name, $start, $end );
                push( @arrayseq, $seq );
                push( @positionss, $genomic_start, $genomic_end );
                $i++;
            }
            close(HITS);
            my @orderbrut1 = sort { $a <=> $b } (@positionss);
            my @order1     = sort { $a <=> $b } (@positionss);
            my $start1     = shift(@order1);
            my $end1       = pop(@order1);
            my $finalseq = join( "NNN", @arrayseq );
            my $FINAL_ID_ID
                = "$header{$query[0]}" . "_" . "$start1" . "_" . "$end1";
            $fasta = Bio::PrimarySeq->new(
                -id  => $FINAL_ID_ID,
                -seq => $finalseq
            );
            $out->write_seq($fasta);
        }

    }

}

#### BED TO FASTA ##

sub BEDTOFASTA {
    my ( $bed, $database ) = @_;
    if ( -e "$database.index" ) { system("rm $database.index"); }
    my $db          = Bio::DB::Fasta->new("$database");
    my $outinterval = Bio::SeqIO->new(
        -file    => ">$bed.fasta",
        - format => "fasta"
    );
    open( INTER, "$bed" );

    while ( my $INTER = <INTER> ) {
        chomp($INTER);
        my @lINTER      = split( " ", $INTER );
        my $CHRin       = $lINTER[0];
        my $sin         = $lINTER[1];
        my $ein         = $lINTER[2];
        my $idin        = join( "_", $CHRin, $sin, $ein );
        my $seqinterval = $db->subseq( $CHRin, $sin, $ein );
        my $fastainter
            = Bio::PrimarySeq->new( -id => $idin, -seq => $seqinterval );
        $outinterval->write_seq($fastainter);
    }
    close(INTER);
}
###### RENAME #####

sub RENAME {
    my ($file) = @_;
    my $seqout = Bio::SeqIO->new(
        -file     => ">$file.rename",
        '-format' => 'Fasta'
    );
    my $i = 1;
    my $seq_in = Bio::SeqIO->new( -format => 'fasta', -file => "$file" );
    while ( my $seq = $seq_in->next_seq() ) {
        my $sequence = $seq->seq;
        my $id       = $seq->primary_id;
        $header{$i} = $id;
        my $fasta = Bio::PrimarySeq->new( -id => $i, -seq => $sequence );
        $seqout->write_seq($fasta);
        $i++;
    }
}
__END__


I) Run TASR 

nohup perl TASR_v.1.0.pl -ref path_to_your_genome/genome.fasta -sfile path_to_your_siRNAs/siRNAs.fasta >& TASR.log &


* to increase speed with parallel threads please use the option -cpu to provide the number of theads used for bowtie2 and BLAST.


II) Example data  

in the folder example_data you will find : 

- genome.fasta (contain Arabidopsis Chr1)
- 24.sirna.fasta (Arabidopsis non redundant 24 nt-siRNAs )

Command line :

nohup perl TASR_v.1.0.pl -ref example_data/genome.fasta -sfile /siRNAs.fasta -usearchv usearch_file_name >& TASR.log &

output file :

the final output will be written to the folder named TASR_output (to change the directory name use -outfold option). This folder will contain multifasta files of different TEs paralogs (by default at least 2 copies, see -cnumber). 


III) the complete list of options is listed bellow : 

* the different possible options are listed bellow (perl TASR_v.1.0.pl)


	TASR version v.1.0

usage: perl TASR_v.1.0.pl -ref genome.fasta -sfile siRNAs.fasta -usearchv path/usearch7 [options]



	-ref          (mandatory)
	 		   	  reference genome sequence in fasta format
		
	-sfile        (mandatory)
			  	  siRNAs sequences in fasta format

        -usearchv     (mandatory)
			 	  please specify the path to executable file name (ex : -usearchv /usr/bin/usearch7.0.1090_i86linux32)
	
	-outfold      (optional)
			 	  output directory (default : -outfold TASR_output)

	-cpu          (optional)
				  minimum number of threads used for mapping and blast (default : -cpu 1)

	-nsirna       (optional)
			 	  minimum number of siRNAs to consider TEs candidates loci (default : -nsirna 4)

	-win          (optional)
				  maximum window size in bp to merge mapped siRNAs (default : -win 150)

	-minlen       (optional)
				  minimum interval length in bp to consider TEs candidates loci default : -minlen 80)
	
	-maxlen       (optional)
				  maximum interval length in bp to consider TEs candidates loci (default : -maxlen 20000)

	-idclust      (optional)
				  minimum % identity to accept blast hits for building families (in [0,1]) (default : -idclust 0.8)

	-overlapclust (optional)
				  minimum % overlap to accept blast hits for building families (in [0,1]) (default : -overlapclust 0.8)

	-cnumber      (optional)
				  minimum copy number to consider TEs families (default : -cnumber 2 )

	-idenblast    (optional)
				  minimum identity percentage cut-off for blast (default : -idenblast 70)

	-evalue       (optional)
				  Expected value (default : -evalue 1e-20)

	-Blastv       (optional)
			 	  If using older blast version please specify '-Blastv NCBI-BLAST' (default : -Blastv NCBI-BLAST+)

	-temp	      (optional)
	 		          temporary files : keep or delete (default : -temp delete)
      
	-h --help
				  print this menu


* For more details about the choice of paramaters please refer the paper.
