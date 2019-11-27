#!/usr/bin/perl
# blast92gff3.pl

=head1 NAME
  
  blast92gff3.pl 
  unpacked that is " BLAST tabular output (-m 9 or 8) conversion to GFF version 3 format "
  
=head1 SYNOPSIS

  ncbi/blastall -m 9 -p tblastn -i daphnia_genes.aa -d aphid-genome -o aphid-daphnia_genes.tblastn
  cat aphid-daphnia_genes.tblastn | sort -k1,1 | blast92gff2.pl > aphid-daphnia_genes.gff

  Note: input must be sorted by query ID.  Default .tblastn result is sorted by query.

=head2 BlastX and tBlastN
  
  This should work the same for BlastX, if you sort on column 2 and use -swap:
     blastall -p blastx -d genes.aa -i genome  | sort -k2,2 | blast92gff2.pl -swap
  It should work for DNA queries like EST/cDNA mapped to a genome also.
  
=head1 ABOUT

  Separate BLAST tabular output (-m 9 or 8) into gene models (match,HSP) by parsing distinct, 
  duplicate gene matches using the query locations as well as source locations.

  This works on one query gene at a time.  The result likely will have many overlapped
  different genes of declining quality.  You can filter those in a 2nd step:
   $td/overbestgene2.perl -in aphid-daphnia_genes.gff > aphid-daphnia_genes-best.gff
 
  Learn more on the subject of accurately locating duplicate genes at 
  http://eugenes.org/gmod/tandy/
    
  See this note on how BioPerl combines distinct gene matches from blat, blast in its
  searchio module
  http://www.bioperl.org/wiki/Talk:GFF_code_audit


=head2 TEST data

  lots of "insect" cuticle tandem duplicate genes in
  daphnia_pulex/scaffold_14:1294000..1322000
  
  http://insects.eugenes.org/cgi-bin/gbrowsenew/gbrowse/daphnia_pulex/?name=scaffold_14:1294000..1322000

  Pick the protein of one cuticle gene,  daphnia:NCBI_GNO_546144
  Run tblastn at wfleabase.org of NCBI_GNO_546144 protein x genome, saving tabular blast result
  blast92gff2.pl < result.tblastn

=head1 AUTHOR

  Don Gilbert, gilbertd@indiana.edu  June 2008
  
=head1 NOTES

=cut


use strict;

use constant MAX_EXON_SEPARATION => 10000; # 50000;  # make an option
use constant QBASEOVER => 5; # query hsp overlap slop in bases
# use constant GFF_LOC => 1;

my $OVERLAP_SLOP_QUERY  = 0.15; # was 0.15; for protein query HSPs
my $OVERLAP_SLOP_GENOME = 0.05; # was 0.50; for genome source HSPs
my $BIN10000= 10000; # binsize
my $GAPSIZE = 400; # what ?? for nearover1() test

my $LOWSCORE_SKIP  = 0.20; # i.e. skip matchs < 50% of max match score

use Getopt::Long;

use vars qw(
$swap_querytarget
$faprefix $debug 
$querySource $exonType
$geneType  $dotarget
$max_bit_score
$bitcutoff $stringent2ndary
%besthash %sumhash @qparthsps @sourcehsps
@mainexons @secondexons @allsaved %moregenes $species $queryid $npart
%shredhash %tandhash $tophsp $lqid $lsid $nwarn $outfile
$gffout %didid
);


my $addgene= 1;
$geneType="match"; # or "protein_match" ...; 
$exonType="HSP"; # or "match_part"; 
$querySource="blast"; 
$faprefix= 'gnl\|'; # drop NCBI extra id prefix
$stringent2ndary=1;
$dotarget= 1;

#  do own sorting for score, loc : still need queryID sort ...
#  ** WARNING : input blast must be sorted by queryID, and best score (sort -k1,1 -k12,12nr) **
#  cat  modelproteins-mygenome.tblastn| sort -k1,1 -k12,12nr |\

my $USAGE =<<"USAGE";
   ** WARNING : input blast must be sorted by queryID  (sort -k1,1 ) **
 cat  modelproteins-mygenome.tblastn| sort -k1,1 |\
      blast92gff3.pl [options]  > modelproteins-mygenome.gff
 options: -swap : swap query, target
 -source $querySource : gff.source field
 -exonType $exonType , -geneType $geneType  : gff feature types
 -lowscore $LOWSCORE_SKIP : skip matches < $LOWSCORE_SKIP of max score
 -notarget, -nomatch  :  change gff output

USAGE

my $ok=&GetOptions( 
  "output=s" => \$outfile,
  "qoverlap=s" => \$OVERLAP_SLOP_QUERY,  
  "overlap=s" => \$OVERLAP_SLOP_GENOME,  
  "LOWSCORE_SKIP=s" => \$LOWSCORE_SKIP,  
  "source=s", \$querySource, # also == GFF source ?
  "geneType=s", \$geneType,  
  "exonType=s", \$exonType,  
  "swap_querytarget!" => \$swap_querytarget,
  'stringent2ndary!' => \$stringent2ndary,
  "faprefix=s", \$faprefix,
  "target!" => \$dotarget,
  "match!" => \$addgene,
  "debug!" => \$debug,
  );
die $USAGE unless($ok);


sub reset_target_vars {
  printGeneLocations(@_);

  $npart= 0; 
  $max_bit_score= 0;
  $bitcutoff=0;
  %shredhash=(); 
  %tandhash=(); #? keep? this ids query regardless of target
  $tophsp='';
  @qparthsps= @sourcehsps= @allsaved= @mainexons= @secondexons= ();
  %moregenes=();
}


my $hascomm=0;
$gffout = *STDOUT;
if( $outfile && open(OUTH,">$outfile")) { $gffout= *OUTH; }
print $gffout "##gff-version 3\n"; ## if($dolocs);

while(<>) {
  chomp;

  if(/^\w/) {  # dont assume blast input has comments .. -m8
    my @v= split "\t";
    my ($qid, $sid, $pctident, $alignment_length, $mismatches, $gap_openings, 
       $q_start, $q_end, $s_start, $s_end, $prob, $bit_score ) = @v;
    next unless($bit_score); # some error logs mixed in.. ## warn, die ??
     
    cleanid($qid);  
    cleanid($sid);  
    if($swap_querytarget) {
      ($qid,$sid,$q_start,$q_end,$s_start,$s_end)= 
      ($sid,$qid,$s_start,$s_end,$q_start,$q_end);
      }

if(1) { # query/gene batching, not target/scaffold, with sort -k1,1 input
    reset_target_vars($lqid) if($qid ne $lqid);
} else {
    reset_target_vars($lsid) if($sid ne $lsid);  
}

    my($s_strand,$q_strand)= ('+','+');
    if ($s_start > $s_end) { $s_strand='-'; ($s_start,$s_end)= ($s_end,$s_start);  }
    if ($q_start > $q_end) { $q_strand='-'; ($q_start,$q_end)= ($q_end,$q_start);  }
    
#     my($sh_id,$sh_start,$sh_end)= deshredloc($sid);
#     $s_start= $sh_start + $s_start;
#     $s_end  = $sh_start + $s_end;
#     $sid    = $sh_id;
    
#     my $skey="$sid.$q_start.$s_start.$s_end.$bit_score";
#     next if($shredhash{$skey});
#     $shredhash{$skey}++;
    
    # if ( $bitcutoff == 0 ) { $bitcutoff= $bit_score * 0.75; } #? 1st == best; also  $npart == 0 
    # elsif ( $stringent2ndary && $bit_score < $bitcutoff && ($sid ne $lsid) ) { next; }
    
    ## binned top hits for filterBesthits ; can we do this inline ?
    # NOT USED NOW# storeBesthits($qid,$sid,$s_start,$s_end,$bit_score);

    # count hsps same query loc, diff db loc
    $tophsp= $sid if ($npart==0);

    my $tkey= "$qid-HSP:$q_start,$q_end";  
    my $saved= 0;
    my $qoverlapped=0;
    my $soverlapped=0;  ## need only for genelocs.other ?

# if(0) { # need this anymore ?
#     foreach my $ex (@qparthsps) {
#       my($bs,$be)= ($$ex[3],$$ex[4]);
#       $qoverlapped= overlap1($q_start, $q_end, $bs, $be, $OVERLAP_SLOP_QUERY);
#       if ($qoverlapped) { $tkey="$qid-HSP:$bs,$be"; last; }
#       }
# 
#     foreach my $ex (@sourcehsps) {
#       unless($$ex[1] ne $sid || $$ex[0] ne $qid) {
#         my($bs,$be)= ($$ex[5],$$ex[6]);
#         $soverlapped= overlap1($s_start, $s_end, $bs, $be, $OVERLAP_SLOP_GENOME);
#         last if $soverlapped;
#         }
#       }
# }

    # this is where we assume input is sorted by query-id and top bitscore : 
    # we drop overlapped hsps silently here ... should reinstitute stringent2ndary test
    next if( $qoverlapped && $soverlapped ); # this happens only for same query, same scaffold
    ##next if($stringent2ndary && $bit_score < $bitcutoff && ($sid ne $lsid));

    my $hspval= [$qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
      $bit_score,  $prob, $tkey, $s_strand,
      $soverlapped, $qoverlapped];
    push(@qparthsps, $hspval) unless($qoverlapped); 
    push(@sourcehsps, $hspval); # unless($qoverlapped && $soverlapped); # use this as recycled hits/scaffold
    
    $npart++; $lqid= $qid; $lsid= $sid;
  } # blast line
} # while(<>)


reset_target_vars($lqid); ## printGeneLocations();
close($gffout) if($gffout);

if(1) { # debug/verbose / stats of parts saved ...
  my($tpart,@tkeys,$v);
  my @sumkeys= sort keys %sumhash;
  warn "# Summary of HSPs saved\n";
  foreach $tpart (@sumkeys) {
    @tkeys= sort keys %{$sumhash{$tpart}};
    foreach (@tkeys) { $v= $sumhash{$tpart}{$_}; warn "# $tpart $_ = $v\n"; }
  }  
}

# end ...................................................



=head2 yet another BLAST assignBestgene algorithm

  algo1 isnt good enough: sort by qgene, source-loc; join hsp by next-nooverlap 
  new algo3: 
      1. input sort by query-gene, top bitscore
      2. save list of top hsp (score), resort hsp by source-location
      3. for top-hsp down, join any next-door left & right same-strand hsp
          mark these as G1..n, save
      4.  continue w/ next-top hsp, eliminating source-overlapped w/ saved G hsp
  
=head2 original algorithm from blast9protstats (messy)

  # step1: -p loc
  gtar -Ozxf $soc/${dp}?prot9.blout*.tgz | \
   perl -n $bg/blast/blast9protstats.pl -p loc > ! ${dp}prot9.blexons &

  # step2: sort hsps with location bins ; should add this to perl
  # ** need this binning to separate, keep best hits/region
  # makeblexonsort.sh
  foreach ble ($dp*.blexons)
  set blo=`echo $ble | sed -e's/blexons/blexonsort/'`
  echo ====== $ble : $blo
  cat $ble | perl -ne\
  '($id)=m/tkey=(\S+)\-HSP:/; ($r)=m/^\w+:(\S+)/; ($tb,$te)=m/tloc=(\d+),(\d+)/;\
  @v=split; $db=$v[1]; $sbin=10000*int($v[3]/10000);  \
  print join("\t",$db,$r,$sbin,$v[5],$v[3],$v[4],$id,$tb,$te),"\n";'  \
  | sort -k1,2 -k3,3n -k4,4nr \
  > $blo
  echo
  end
  
  # step3: -p table (default)
  gtar -Ozxf $soc/${dp}?prot9.blout*.tgz | \
   perl -n $bg/blast/blast9protstats.pl -over 0.5 -besthits ${dp}prot9.blexonsort  -out ${dp}prot9.bltab7 &

=cut

sub _sortHsp_Score {
  my($ap)= $a->[7];
  my($bp)= $b->[7];
  return ($bp <=> $ap); # top score 1st
}

sub _sortHsp_SrcLocation  {
  my($ar,$ab,$ae,$ao)= @{$a}[1,5,6,10]; 
  my($br,$bb,$be,$bo)= @{$b}[1,5,6,10];
  my $ocmp= ($ao eq $bo)? 0 : ($ao eq "-") ? 1 : -1;
  return ($ar cmp $br || $ocmp || $ab <=> $bb || $be <=> $ae); # sort all -strand  last
}

sub bestlocHspset {
  my($hsploc, $tsid, $ts_start, $ts_end, $ts_strand)= @_;
  my @before=(); my @after=();
  my $skiphsp= 0;
  my($trange0, $trange1)= ($ts_start - MAX_EXON_SEPARATION, $ts_end + MAX_EXON_SEPARATION);
  foreach my $hspval (@$hsploc) { # instead of foreach can we hash-find nearby hsps?
    my($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score,  $prob, $tkey, $s_strand,
        $soverlapped, $qoverlapped)
        = @$hspval;
    next unless($sid eq $tsid && $s_strand eq $ts_strand && $s_start > $trange0 && $s_end < $trange1);
    if($s_end <= $ts_start + QBASEOVER) { 
      unshift(@before, $hspval); 
      }
    elsif($s_start >= $ts_end - QBASEOVER) { 
      push(@after, $hspval); 
      }
    }  
  return (\@before, \@after); # sorted around hspbest
}


sub assignBestgene {  # version 3
  #? my($theqid) = @_; # not used

  my($topsid);
  $npart=0;
  my $saved=0;
  my $genenum = 0; my $lastgenenum= 0;  
  my $lastexon= undef;
  @allsaved=(); # dont need global
  
  my @hspbest = sort _sortHsp_Score @sourcehsps;
  my @hsploc  = sort _sortHsp_SrcLocation @sourcehsps;
  
  # input sourcehsp should all be for same query-gene, sorted by location (NOT/was bitscore)
  # need location-sort to match up exon parts of genes
  my %donehsp=();
  my($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score,  $prob, $tkey, $s_strand,
        $soverlapped, $qoverlapped);
  
foreach my $hspbest (@hspbest) {
  ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score,  $prob, $tkey, $s_strand,
        $soverlapped, $qoverlapped)
        = @$hspbest;
  my $topkey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
# warn "#top $topkey\n"; #DEBUG
  next if $donehsp{$topkey}++;
  next if overlapsome($s_start, $s_end, \@allsaved, $OVERLAP_SLOP_GENOME);
    
  $genenum++;
  $topsid= $sid; ## if ($npart==0);

  my $keynum= $genenum;
  unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
  push( @{$moregenes{$keynum}}, $hspbest); $saved=1;
  $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
  if ($saved) {
    push(@allsaved, $hspbest);
    $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
    }
  
  my($tsid, $tq_start, $tq_end, $ts_start, $ts_end, $ts_strand)= 
    ($sid, $q_start, $q_end, $s_start, $s_end, $s_strand);
  my($trange0, $trange1)= ($ts_start - MAX_EXON_SEPARATION, $ts_end + MAX_EXON_SEPARATION);
  my($srange0, $srange1)= ($ts_start, $ts_end);
  my($qrange0, $qrange1)= ($tq_start, $tq_end);
  
    ## FIXME.3 this really needs to step thru @hsploc starting at $hspbest loc and go down,up from there
    ## otherwise qrange, srange are bad.  are getting two genes made from very good pieces of 1 gene match
    ## due to interior hsp's skipped in first pass nearest to hspbest.
    
  my(@before,@after);
  my($before, $after)= bestlocHspset(\@hsploc, $tsid, $ts_start, $ts_end, $ts_strand);
  
  foreach my $hspval (@$after, @$before) { # instead of foreach can we hash-find nearby hsps?
    # next unless (ref $hspval); #?
    
    ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score,  $prob, $tkey, $s_strand,
        $soverlapped, $qoverlapped)
        = @$hspval;

    ## FIXME here; should not skip, but keep some of these to check; 
    ## should not look far away if nearby exon fits; it it is done already or overlaps, count in qrange and skip on
    ## FIXME.2 new problem with this; far-hsp already done can eat away query-range; need to skip those

    ## already done# next unless($sid eq $tsid && $s_strand eq $ts_strand && $s_start > $trange0 && $s_end < $trange1);

    my $skiphsp= 0;
    my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
## warn "#at $atlockey\n"; #DEBUG
    next if($atlockey eq $topkey);
    $skiphsp=1 if $donehsp{$atlockey};
    $skiphsp=1 if $skiphsp || overlapsome($s_start, $s_end, \@allsaved, $OVERLAP_SLOP_GENOME);

    my $qover= overlap1($q_start, $q_end, $qrange0, $qrange1, $OVERLAP_SLOP_QUERY);
    next if($qover); # last?
    # $skiphsp=1 if ($qover);
    ## FIXME: need to look at qloc vs top-qloc; if -strand, @before must be higher qloc, @after lower
    ## and v.v.

    # if($skiphsp) now check that s_start,s_end is *near* top hsp; skip if not
    my $nearover= nearover1($s_start, $s_end, $srange0, $srange1, $OVERLAP_SLOP_GENOME, 5000); # GAPSIZE
    if($nearover >= 0) { next if $skiphsp; }
    
    if($s_end <= $ts_start + QBASEOVER) { # before
#       if($s_strand eq '-') { $qrange1 = $q_end if($q_end> $qrange1);  }
#       else {  $qrange0 = $q_start if($q_start< $qrange0); } # not before      
      if($s_strand eq '-') { next if($q_start + QBASEOVER < $qrange1); $qrange1 = $q_end if($q_end> $qrange1);  }
      else { next if($q_end - QBASEOVER > $qrange0);  $qrange0 = $q_start if($q_start< $qrange0); } # not before      
      unshift(@before, $hspval) unless $skiphsp; 
      $srange0 = $s_start unless $skiphsp;
      }
      
    elsif($s_start >= $ts_end - QBASEOVER) { # after; bug in next here
#       if($s_strand eq '-') {  $qrange0 = $q_start if($q_start < $qrange0); }
#       else {  $qrange1 = $q_end if($q_end > $qrange1);  } # not after
      if($s_strand eq '-') { next if($q_end - QBASEOVER > $qrange0); $qrange0 = $q_start if($q_start< $qrange0); }
      else { next if($q_start  + QBASEOVER < $qrange1);  $qrange1 = $q_end if($q_end> $qrange1);  } # not after
      push(@after, $hspval) unless $skiphsp; 
      $srange1 = $s_end unless $skiphsp;
      }
##warn "#skip=$skiphsp $atlockey\n"; #DEBUG
  }  
  
  #? limit before, after size? dont try to find what isn't there by too far a match
  foreach my $hspval (@before, @after) { # .. and after
  
    ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
        $bit_score,  $prob, $tkey, $s_strand,
        $soverlapped, $qoverlapped)
        = @$hspval;
    my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
    $saved= 0; ##no## $genenum= 0;
       
    if(1) {
      #my $keynum= "G" . $genenum;
      my $keynum=  $genenum;
      unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
      push( @{$moregenes{$keynum}}, $hspval); $saved=1;
      $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
    }

    if ($saved) {
      $lastexon= $hspval;
      $lastgenenum= $genenum;
      $donehsp{$atlockey}++;
      push(@allsaved, $hspval);
      $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
      }
   $npart++;
  }
}

}

sub printGeneLocations {
  my($qid) = @_; # not used
  
  if(@sourcehsps && $qid && 0 < $didid{$qid}++) { die "ERROR: $qid already seen\n$USAGE"; }
  
  assignBestgene(); # NEW
  
  my $ng= printOneLocation(\@mainexons, 'G1'); # not used now
  $ng  += printOneLocation(\@secondexons, 'G2'); # not used now
  foreach my $keynum (sort{$a<=>$b} keys %moregenes) { ## NOW NUMERIC KEY
    (my $gnum= $keynum) =~ s/^(.)[a-z]+/$1/;  ## S2..S9; o2..o9
    $gnum = 'G'.$gnum;
    $ng += printOneLocation($moregenes{$keynum}, $gnum);
    }
  return $ng;
}


sub printOneLocation {
  my($exons, $igene) = @_;
  return 0 unless(ref $exons && @$exons > 0);

  my $qdb=".";
  my @locs= ();
  my($sum_bit_score,$sum_align,$iex)=(0)x10;
  my($qid1,$qid,$sid,$alignment_length,$q_start,$q_end,$s_start,$s_end,
      $bit_score, $prob, $tkey,$s_strand, $soverlapped, $qoverlapped);

  my @gffout= ();
  my($g_start,$g_end, $g_strand)=(undef,0);
  foreach my $ex ((sort{ $$a[5] <=> $$b[5] } @$exons)) {
    ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
      $bit_score,  $prob, $tkey, $s_strand,
      $soverlapped, $qoverlapped)= @$ex;
    
    $qdb= $querySource; # default
    ($qdb,$qid1)= ($qid =~m/:/) ? split(/[:]/, $qid,2) : ($qdb,$qid);
    
    $g_start= $s_start unless(defined $g_start);
    $g_end= $s_end;
    $g_strand= $s_strand;
    $iex++;
    $sum_bit_score += $bit_score;
    $sum_align += $alignment_length;
    if(1) { # == GFF_LOC
      $tkey =~ s/^\w+://;  $tkey =~ s/,/-/g; # extra dbid and no commas in values
      #was# $attr= "Parent=${qid1}_${igene};tkey=$tkey;tloc=$q_start-$q_end;align=$alignment_length";
      my $attr= "Parent=${qid1}_${igene}";
      $attr.= ";Target=$qid $q_start $q_end;align=$alignment_length" if($dotarget);
      # print $gffout 
      push @gffout, join("\t",
      ($sid, $qdb, $exonType, $s_start, $s_end, $bit_score, $s_strand,".",$attr)). "\n";
      }
    # else { push(@locs,"$s_start..$s_end"); }
    }
    
  $max_bit_score = $sum_bit_score if($sum_bit_score > $max_bit_score ); # assumes output by best match 1st
  return if ($sum_bit_score / $max_bit_score < $LOWSCORE_SKIP);
  
  if($addgene) {
    my $attr= "ID=${qid1}_${igene}";
    $attr.= ";Target=$qid;align=$sum_align" if($dotarget);
    unshift @gffout, join("\t",
      ($sid, $qdb, $geneType, $g_start, $g_end, $sum_bit_score, $g_strand,".",$attr)). "\n";
  }
  
  print $gffout @gffout;
  return 1;
}

sub overlapsome { # need scaffold?
  my($s_start, $s_end, $exons, $OVERLAP_SLOP)= @_; 
  foreach my $ex (@$exons) {
    return 1 if overlap1($s_start, $s_end, $$ex[5], $$ex[6], $OVERLAP_SLOP);
  }
  return 0;
}

sub overlap1 {
  my($s_start, $s_end, $b_start, $b_end, $OVERLAP_SLOP)= @_; 
  return 0 unless ( $s_start < $b_end && $s_end > $b_start  ); ## no overlap
  return 1 if ( $s_start >= $b_start && $s_end <= $b_end ); # contained-in
  if ( $s_end > $b_start && $s_start < $b_end ) {  
    ## e.g.  s= 20,50  ; b = 10,30 : o= 10* or 40
    my $b_len= 1 + $b_end - $b_start;
    my $s_len= 1 + $s_end - $s_start; # choose ? biggest
    my $olp1 =  $b_end - $s_start; # 30-20 = 10
    my $olp2 =  $s_end - $b_start; # 50-10 = 40
    $olp1= $olp2 if ($olp2<$olp1);
    $s_len= $b_len if ($s_len<$b_len);
    return 1 if (($olp1 / $s_len) > $OVERLAP_SLOP);
    }
  return 0;
}  

sub nearover1 {  # near == -1; over == +1; other == 0
  my($s_start, $s_end, $b_start, $b_end, $OVERLAP_SLOP, $gapsize)= @_; 
  
  # need some overlap slop here; note locs here are QUERY == protein; should not be large
  unless ( $s_start < ($b_end - QBASEOVER) && $s_end > ($b_start + QBASEOVER)  ) {
    if($s_start >= $b_end) { return ($s_start - $b_end < $gapsize) ? -1 : 0; }
    elsif($b_start >= $s_end) { return ($b_start - $s_end < $gapsize) ? -1 : 0; }
    return 0;
  }
  
  return overlap1($s_start, $s_end, $b_start, $b_end, $OVERLAP_SLOP);
}  


sub minabs {
  my($a,$b)= @_;
  $a= abs($a); $b= abs($b);
  return ($a>$b)? $b : $a;
}

sub toofar1 {
  my($q_start, $q_end, $b_start, $b_end)= @_; 
  return 1 if minabs($b_end - $q_start,$q_end - $b_start) > MAX_EXON_SEPARATION;
  return 0;
}

sub cleanid { 
  unless( $_[0] =~ s/$faprefix//) { $_[0] =~ s/^gi\|\d+\|(\S+)/$1/; }  # drop gi nums for ids
  $_[0] =~ s/\|$//; ## some ids have '|' at end of id; chomp
  $_[0] =~ s/\|/:/; #?? change pipe to db:id format
}

__END__

my $bestsort=0; # fixme

# NOT USED NOW
sub storeBesthits {
  my($qid,$sid,$s_start,$s_end,$bit_score)= @_;
  
  ## binned top hits for filterBesthits ; can we do this inline ?
  my($qdb,$qid1)= ($qid =~m/:/) ? split(/[:]/, $qid,2) : ($querySource,$qid);
  my ($sbin)= ($BIN10000 * int($s_start/$BIN10000)) ; # must match besthash bin size
  my $bestkey="$qdb.$sid.$sbin";
  my @bv= ($bit_score,$s_start,$s_end,$qid);  
  push( @{$besthash{$bestkey}}, \@bv); # arrayref of arrayref
  
  $bestkey="$qdb"; # best over all genome == _G1
  if( ! $besthash{$bestkey} ) {
    $besthash{$bestkey}= [\@bv];
  } elsif($bit_score > $besthash{$bestkey}->[0]->[0]) {
    $besthash{$bestkey}->[0]= \@bv;
  }
}

# NOT USED NOW
sub filterBesthits {
  my($qid,$sid,$s_start,$s_end,$s_strand,$bit_score)= @_;
 
  unless($bestsort>0) {
    $bestsort++;
    foreach my $akey (sort keys %besthash) {
      my @blocs= sort{$b->[0] <=> $a->[0]} @ { $besthash{$akey} };
      $besthash{$akey}= \@blocs;
      }
    }
    
  my $snotbest= 0;
  return 0 unless(scalar(%besthash));
 
  my ($qdb,$qid1)= ($qid =~m/:/) ? split(/[:]/, $qid,2) : ($querySource,$qid);
  my ($sbin)= ($BIN10000 * int($s_start/$BIN10000)) ; # must match besthash bin size
  my($sdb,$sid1)= ($sid =~/:/) ? split(/[:]/, $sid,2) : ('',$sid);
  my $bestkey="$qdb.$sid1.$sbin";
  my $isbest = $besthash{$bestkey};

  if($isbest) {
    my @blocs= @ { $besthash{$bestkey} };
    foreach my $bl (@blocs) {
      my($bbits,$bb,$be,$bqid)= @$bl;  
      $snotbest= overlap1($s_start, $s_end, $bb, $be, $OVERLAP_SLOP_GENOME);
      if($snotbest) {
        $snotbest= 0 if($bit_score >= $bbits || $qid eq $bqid); # keep self match
        last;
        }
      }
    $sumhash{'besthits'}{($snotbest?'no':'yes')}++; 
    }
  else {
    $sumhash{'besthits'}{'missed'}++; 
    }
  return $snotbest;
}


# NOT USED NOW
sub pushGeneExons {
  my($exons, $hspval, $q_start,$q_end,$s_start,$s_end) = @_;

  my $qoverlapped=0;
  my $toofar= -1;
  foreach my $ex (@$exons) {
    my ($bs,$be)= ($$ex[3],$$ex[4]);
    $qoverlapped= overlap1($q_start,$q_end, $bs, $be, $OVERLAP_SLOP_QUERY);
    last if ($qoverlapped);
    if($toofar != 0) { $toofar= toofar1($s_start,$s_end, $$ex[5],$$ex[6]); }
    }
  return 0 if($qoverlapped || $toofar>0);
  push(@$exons, $hspval); # [$qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, $bit_score,  $prob]
  return 1; 
}




## this one is out; not good enough
# 
# sub assignBestgene1 {
#   my($tophsp);
#   $npart=0;
#   my $saved=0;
#   my $genenum = 0; my $lastgenenum= 0;  
#   my $lastexon= undef;
# 
#   # input sourcehsp should all be for same query-gene, sorted by location (NOT/was bitscore)
#   # need location-sort to match up exon parts of genes
#   
#   foreach my $hspval (@sourcehsps) { 
#     my($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped)
#         = @$hspval;
#     $saved= 0; $genenum= 0;
#      
# ## test this; not working right
# ## really need to test, for each hsp?, nearby hsp both directions and need to keep strand same
#     my $overlast=0; 
#     foreach my $ex ($lastexon) {
#       my ($bsid,$bs,$be)= ($$ex[1],$$ex[3],$$ex[4]);
#       next unless($bsid eq $sid);
#       # $overlast= overlap1($q_start, $q_end, $bs, $be, $OVERLAP_SLOP_QUERY);
#       # last if ($overlast);
#       my $nearlast= nearover1($q_start, $q_end, $bs, $be, $OVERLAP_SLOP_QUERY);
#       if ($nearlast<0) { $genenum= $lastgenenum; last; }
#       elsif($nearlast>0) { $overlast=1; last; }
#       }
#     # next if $overlast; #??
#   
#     $tophsp= $sid if ($npart==0);
#     my $issame= ($sid eq $tophsp); 
#     my $scafpart;
#     if($npart==0 || $issame) { $tandhash{'Samescaf'}{$tkey}++; $scafpart='Samescaf'; }
#     if($npart==0 || !$issame) { $tandhash{'otherscaf'}{$tkey}++; $scafpart= 'otherscaf' if($npart>0); }
# 
#     $genenum = $tandhash{$scafpart}{$tkey} unless ($genenum); 
#         # ^^ this is easily bad : tkey and messes up joining exon-parts of gene : splits too much
#     
#     if( $soverlapped  ) { 
#       $sumhash{'other_soverlap'}{($saved?'saved':'notsaved')}++; 
# 
#     #} elsif( $overlast  ) { # this is problem
#     #  $sumhash{'other_qoverlap'}{($saved?'saved':'notsaved')}++; 
#     
#     } elsif($sid eq $tophsp && $genenum == 1) {  #  primary-match exons
#       my $snotbest= 0; #?? filterBesthits($qid,$sid,$s_start,$s_end,$s_strand,$bit_score); 
#       if ($snotbest) { $saved=0; }
#       else { $saved= pushGeneExons(\@mainexons, $hspval, $q_start,$q_end, $s_start,$s_end); }
#       $sumhash{'G1'}{($saved?'saved':'notsaved')}++; 
#       }
# 
# ## test off; need this?      
# #     elsif($sid eq $tophsp && $genenum == 2) { #  2nd-match exons
# #       my $snotbest= filterBesthits($qid,$sid,$s_start,$s_end,$s_strand,$bit_score);  
# #       if ($snotbest) { $saved=0; }
# #       else { $saved= pushGeneExons(\@secondexons, $hspval, $q_start,$q_end, $s_start,$s_end); }
# #       $sumhash{'G2'}{($saved?'saved':'notsaved')}++; 
# #       }
# 
#     elsif( !$soverlapped ) { # other-match exons
#       my $snotbest= filterBesthits($qid,$sid,$s_start,$s_end,$s_strand,$bit_score);  
#       $saved= 0;
#       unless($snotbest) {
#         my $keynum= $scafpart . $genenum;
#         unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
#         $saved= pushGeneExons( $moregenes{$keynum}, $hspval, $q_start,$q_end, $s_start,$s_end);
#       }
#       $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
#     }
# 
#     if ($saved) {
#       $lastexon= $hspval;
#       $lastgenenum= $genenum;
#       push(@allsaved, $hspval);
#       $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
#       }
#    $npart++;
#   }
# }


## this one is out; not good enough
# sub assignBestgene_OLD {
#   my($tophsp);
#   $npart=0;
#   my $saved=0;
#   foreach my $hspval (@sourcehsps) { # these should all be for same query-gene, sorted by bitscore
#      my($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped)
#         = @$hspval;
#   
#     $tophsp= $sid if ($npart==0);
#     my $issame= ($sid eq $tophsp); 
#     my $scafpart;
#     if($npart==0 || $issame) { $tandhash{'Samescaf'}{$tkey}++; $scafpart='Samescaf'; }
#     if($npart==0 || !$issame) { $tandhash{'otherscaf'}{$tkey}++; $scafpart= 'otherscaf' if($npart>0); }
# 
#     my $genenum = $tandhash{$scafpart}{$tkey}; ## need same/diff num
#         # ^^ this is easily bad : tkey and messes up joining exon-parts of gene : splits too much
#     
#     if($sid eq $tophsp && $genenum == 1) {  #  primary-match exons
#       my $snotbest= 0; #?? filterBesthits($qid,$sid,$s_start,$s_end,$s_strand,$bit_score); 
#       if ($snotbest) { $saved=0; }
#       else { $saved= pushGeneExons(\@mainexons, $hspval, $q_start,$q_end, $s_start,$s_end); }
#       $sumhash{'G1'}{($saved?'saved':'notsaved')}++; 
#       }
#       
#     elsif($sid eq $tophsp && $genenum == 2) { #  2nd-match exons
#       my $snotbest= filterBesthits($qid,$sid,$s_start,$s_end,$s_strand,$bit_score);  
#       if ($snotbest) { $saved=0; }
#       else { $saved= pushGeneExons(\@secondexons, $hspval, $q_start,$q_end, $s_start,$s_end); }
#       $sumhash{'G2'}{($saved?'saved':'notsaved')}++; 
#       }
# 
#    elsif( !$soverlapped) { # other-match exons
#       my $snotbest= filterBesthits($qid,$sid,$s_start,$s_end,$s_strand,$bit_score);  
#       $saved= 0;
#       unless($snotbest) {
#       my $keynum= $scafpart . $genenum;
#       unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
#       $saved= pushGeneExons( $moregenes{$keynum}, $hspval, $q_start,$q_end, $s_start,$s_end);
#       }
#       $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
#     }
# 
#     if ($saved) {
#       push(@allsaved, $hspval);
#       $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
#       }
#    $npart++;
#   }
# }
# 
    

# sub assignBestgene2 {  # version 2
# 
#   my($topsid);
#   $npart=0;
#   my $saved=0;
#   my $genenum = 0; my $lastgenenum= 0;  
#   my $lastexon= undef;
#   @allsaved=(); # dont need global
#   
#   my @hspbest = sort _sortHsp_Score @sourcehsps;
#   my @hsploc  = sort _sortHsp_SrcLocation @sourcehsps;
#   
#   # input sourcehsp should all be for same query-gene, sorted by location (NOT/was bitscore)
#   # need location-sort to match up exon parts of genes
#   my %donehsp=();
#   my($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped);
#   
# foreach my $hspbest (@hspbest) {
#   ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped)
#         = @$hspbest;
#   my $topkey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
#   next if $donehsp{$topkey}++;
#   next if overlapsome($s_start, $s_end, \@allsaved, $OVERLAP_SLOP_GENOME);
#     
#   $genenum++;
#   $topsid= $sid; ## if ($npart==0);
# 
#   my $keynum= $genenum;
#   unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
#   push( @{$moregenes{$keynum}}, $hspbest); $saved=1;
#   $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
#   if ($saved) {
#     push(@allsaved, $hspbest);
#     $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
#     }
#   
#   my($tsid, $tq_start, $tq_end, $ts_start, $ts_end, $ts_strand)= 
#     ($sid, $q_start, $q_end, $s_start, $s_end, $s_strand);
#   my($trange0, $trange1)= ($ts_start - MAX_EXON_SEPARATION, $ts_end + MAX_EXON_SEPARATION);
#   my($srange0, $srange1)= ($ts_start, $ts_end);
#   my($qrange0, $qrange1)= ($tq_start, $tq_end);
#   my(@before,@after);
#   
#   foreach my $hspval (@hsploc) { # instead of foreach can we hash-find nearby hsps?
#   
#     ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped)
#         = @$hspval;
# 
#     ## FIXME here; should not skip, but keep some of these to check; 
#     ## should not look far away if nearby exon fits; it it is done already or overlaps, count in qrange and skip on
#     ## FIXME.2 new problem with this; far-hsp already done can eat away query-range; need to skip those
# 
#     ## FIXME.3 this really needs to step thru @hsploc starting at $hspbest loc and go down,up from there
#     ## otherwise qrange, srange are bad.  are getting two genes made from very good pieces of 1 gene match
#     ## due to interior hsp's skipped in first pass nearest to hspbest.
#     
#     
#     my $skiphsp= 0;
#     next unless($sid eq $tsid && $s_strand eq $ts_strand && $s_start > $trange0 && $s_end < $trange1);
# 
#     my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
#     $skiphsp=1 if $donehsp{$atlockey};
#     $skiphsp=1 if $skiphsp || overlapsome($s_start, $s_end, \@allsaved, $OVERLAP_SLOP_GENOME);
# 
#     my $qover= overlap1($q_start, $q_end, $qrange0, $qrange1, $OVERLAP_SLOP_QUERY);
#     next if($qover); # last?
#     # $skiphsp=1 if ($qover);
#     ## FIXME: need to look at qloc vs top-qloc; if -strand, @before must be higher qloc, @after lower
#     ## and v.v.
# 
#     # if($skiphsp) now check that s_start,s_end is *near* top hsp; skip if not
#     my $nearover= nearover1($s_start, $s_end, $srange0, $srange1, $OVERLAP_SLOP_GENOME, 1000); # GAPSIZE
#     if($nearover >= 0) { next if $skiphsp; }
#     
#     if($s_end <= $ts_start + QBASEOVER) { 
#       if($s_strand eq '-') { next if($q_end < $qrange1); $qrange1 = $q_end if($q_end> $qrange1);  }
#       else { next if($q_start > $qrange0);  $qrange0 = $q_start if($q_start< $qrange0); } # not before      
#       unshift(@before, $hspval) unless $skiphsp; 
#       $srange0 = $s_start unless $skiphsp;
#       }
#     elsif($s_start >= $ts_end - QBASEOVER) { 
#       if($s_strand eq '-') { next if($q_start > $qrange0); $qrange0 = $q_start if($q_start< $qrange0); }
#       else { next if($q_end < $qrange1);  $qrange1 = $q_end if($q_end> $qrange1);  } # not after
#       push(@after, $hspval) unless $skiphsp; 
#       $srange1 = $s_end unless $skiphsp;
#       }
#   }  
#   
#   #? limit before, after size? dont try to find what isn't there by too far a match
#   foreach my $hspval (@before, @after) { # .. and after
#   
#     ($qid,$sid,$alignment_length, $q_start,$q_end,$s_start,$s_end, 
#         $bit_score,  $prob, $tkey, $s_strand,
#         $soverlapped, $qoverlapped)
#         = @$hspval;
#     my $atlockey="$qid.$sid.$q_start.$s_start.$s_end.$bit_score";
#     $saved= 0; ##no## $genenum= 0;
#        
#     if(1) {
#       #my $keynum= "G" . $genenum;
#       my $keynum=  $genenum;
#       unless(ref $moregenes{$keynum}) { $moregenes{$keynum}= []; }
#       push( @{$moregenes{$keynum}}, $hspval); $saved=1;
#       $sumhash{'other'}{($saved?'saved':'notsaved')}++; 
#     }
# 
#     if ($saved) {
#       $lastexon= $hspval;
#       $lastgenenum= $genenum;
#       $donehsp{$atlockey}++;
#       push(@allsaved, $hspval);
#       $sumhash{'ALL'}{($saved?'saved':'notsaved')}++; 
#       }
#    $npart++;
#   }
# }
# 
# }


