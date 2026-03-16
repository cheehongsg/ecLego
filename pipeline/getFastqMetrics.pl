#!/bin/env perl

##
## ecLego3 v3.2403.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##
## ecLego component
##   get sequencing metrics and read length distribution
##   useful for assess suitable parameters for assembly and polishing
##

use strict;
use Data::Dumper;

my $fqFile = $ARGV[0];
my $minLength = 0;
$minLength = int($ARGV[1]) if (defined $ARGV[1]);
$minLength = 0 if ($minLength<0);

if ($fqFile =~ /\.gz$/) {
  open INFILE, "gzip -dc $fqFile |" || die "Fail to open $fqFile\n$!\n";
} else {
  open INFILE, $fqFile || die "Fail to open $fqFile\n$!\n";
}

my @readLens = ();
my $totalLength = 0;
my $numTotalReads = 0;
my $numPassReads = 0;
my $totalPassLength = 0;
while (my $id=<INFILE>) {
  my $seq = <INFILE>; chomp();
  my $seqLen = length($seq);
  if ($seqLen>$minLength) {
    $numPassReads++;
    $totalPassLength += $seqLen;
  }
  $numTotalReads++;
  $totalLength += $seqLen;
  push @readLens, $seqLen;

  my $ignore = <INFILE>;
  $ignore = <INFILE>;
}
close INFILE;

printf "# Total read length       : %d\n", $totalLength;
printf "# Total reads             : %d\n", $numTotalReads;
if ($minLength>0) {
  printf "# Minimum length threshold: %d\n", $minLength;
  printf "# Passed read length      : %d\n", $totalPassLength;
  printf "# Passed reads            : %d\n", $numPassReads;
}

# calculate a series of nX: N10, N30, N50, N70, N90
@readLens = sort { $b<=>$a } @readLens;
my @Nxs = ();
push @Nxs, { id=>'N10', lenThreshold=>(0.1 * $totalLength), numReads=>0, lenRead=>0};
push @Nxs, { id=>'N30', lenThreshold=>(0.3 * $totalLength), numReads=>0, lenRead=>0};
push @Nxs, { id=>'N50', lenThreshold=>(0.5 * $totalLength), numReads=>0, lenRead=>0};
push @Nxs, { id=>'N70', lenThreshold=>(0.7 * $totalLength), numReads=>0, lenRead=>0};
push @Nxs, { id=>'N90', lenThreshold=>(0.9 * $totalLength), numReads=>0, lenRead=>0};
@Nxs = reverse @Nxs;

@readLens = sort { $b<=>$a } @readLens;
my $numReads = 0;
my $cummulativeLen = 0;
foreach my $seqLen (@readLens) {
  $numReads++;
  $cummulativeLen += $seqLen;
  foreach my $nxsRef (@Nxs) {
    if (0==$nxsRef->{numReads} && $cummulativeLen > $nxsRef->{lenThreshold}) {
      $nxsRef->{numReads} = $numReads;
      $nxsRef->{lenRead} = $seqLen;
    }
  }
}

my $n90Ref = undef;
foreach my $nxsRef (reverse @Nxs) {
  printf "# %s = %d (%d reads, %.2f%%)\n", 
    $nxsRef->{id}, $nxsRef->{lenRead}, $nxsRef->{numReads},
    $nxsRef->{numReads} * 100.0 / $numTotalReads;
  if ('N90' eq $nxsRef->{id}) {
    $n90Ref = $nxsRef;
  }
}
if (defined $n90Ref) {
  if ($n90Ref->{lenRead} < 10000) {
    if ($n90Ref->{lenRead} < 1000) {
      print "# ERROR: N90 is less than 1 kb, which make assembly impossible.\n";
    } elsif ($n90Ref->{lenRead} < 3500) {
      print "# WARNING: N90 is less than 3.5 kb, which make assembly very difficult.\n";
      my $recommendation = int($n90Ref->{lenRead} / 500) * 500;
      print "# RECOMMENDATION: Consider using a minimum overlap threshold of $recommendation.\n";
    } else {
      print "# WARNING: N90 is less than 10kb, which will impact assembly.\n";
      my $recommendation = int($n90Ref->{lenRead} / 500) * 500;
      print "# RECOMMENDATION: Consider using a minimum overlap threshold of $recommendation.\n";
    }
  }
}


# TODO: parameterize the effective genome size
my $effectiveGenomeSize = 3.2 * 1e9;
print "#\n";
printf "# Estimated coverage      : %.2fx\n", $totalLength/$effectiveGenomeSize;
printf "# Please set the number of reads per two-copy genome to %d\n", int($totalLength/$effectiveGenomeSize);

exit 0;
