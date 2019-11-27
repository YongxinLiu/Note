#!/usr/bin/perl
use warnings;
use strict;

while (<>) {

  unless (/^\@/) {
    warn "$_ should have had an \@ at the start and it didn't\n";
    next;
  }
  my $id1 = $_;
  my $seq = <>;
  my $id2 = <>;
  my $qual = <>;

  if ($seq =~/^[@+]/) {
    warn "Sequence '$seq' looked like an id";
    next;
  }
  if ($qual =~/^[@+]/) {
    warn "Quality '$qual' looked like an id";
    next;
  }
  if ($id2 !~ /^\+/) {
    warn "Midline '$id2' didn't start with a +";
    next;
  }

  if ($qual =~ /[GATCN]{20,}/) {
    warn "Quality '$qual' looked like sequence";
    next;
  }

  if (length($seq) != length($qual)) {
    warn "Seq $seq and Qual $qual weren't the same length";
    next;
  }

  print $id1,$seq,$id2,$qual;


}