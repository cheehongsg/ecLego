#!/bin/env perl

##
## ecLego3 v3.2304.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##
## ecLego component
##   calculate coverage thresholds based on k-mer frequency distribution
##

use strict;
use Data::Dumper;
use POSIX qw/floor ceil/;
my $G_WindowFlank = 2;

my $kmerGraphFile = $ARGV[0];
open INFILE, $kmerGraphFile || die "Fail to open $kmerGraphFile\n$!\n";
my @counts = ();
while (<INFILE>) {
  chomp();
  my ($kmerCount, $kmerFrequency) = split(/\t/);
  push @counts, {kmercount=>$kmerCount, frequency=>$kmerFrequency};
 
}
close INFILE;

my $cov2CN = floor(get2CNCoverage(\@counts));
my $cov1CN = floor($cov2CN/2);
my $cov4CN = floor($cov2CN*2);
my $cov5CN = floor($cov2CN*5/2);
my $cov6CN = floor($cov2CN*3);
my $errorRate = estimateErrorRate(\@counts);

printf "# %s\n", $kmerGraphFile;
printf "covCN=%d\n", $cov1CN;
printf "cov2CN=%d\n", $cov2CN;
printf "cov4CN=%d\n", $cov4CN;
printf "cov5CN=%d\n", $cov5CN;
printf "cov6CN=%d\n", $cov6CN;
printf "errorRate=%.3f\n", $errorRate;

exit 0;

sub getFirstValleyIndex {
  my ($dataPointsRef) = @_;
  my $numDatapoints = scalar(@{$dataPointsRef});
  my $valleyIndex = 0;
  my $valleyRef = $dataPointsRef->[$valleyIndex];
  for(my $i=1; $i<($numDatapoints-1); ++$i) {
    if ($valleyRef->{frequency}>$dataPointsRef->[$i]->{frequency}) {
      $valleyRef = $dataPointsRef->[$i];
      $valleyIndex = $i;
    } else {
      last;
    }
  }
  return $valleyIndex;
}

sub get2CNCoverage {
  my ($dataPointsRef) = @_;

  # find first valley
  my $startIndex = getFirstValleyIndex($dataPointsRef);

  # window averaging
  my $numDatapoints = scalar(@{$dataPointsRef});
  my $windowFlank = $G_WindowFlank;
  my $numElementsInWindow = 1+2*$windowFlank;
  my @winAves = ();
  for(my $i=$startIndex+$windowFlank; $i<($numDatapoints-$windowFlank); ++$i) {
    my $kmerCount = $dataPointsRef->[$i]->{kmercount};
    my $totalFreq = 0;
    for(my $j=$i-$windowFlank; $j<=($i+$windowFlank); ++$j) {
        $totalFreq += $dataPointsRef->[$i-1]->{frequency};
    }
    push @winAves, {kmercount=>$kmerCount, frequency=>int($totalFreq/$numElementsInWindow)};
  }

  # find first peak
  my $numWindowPoints = scalar(@winAves);
  for(my $i=$windowFlank; $i<($numWindowPoints-$windowFlank); ++$i) {
    my $failed = 0;
    for(my $j=$i; $j<($i+$windowFlank); ++$j) {
      if ($winAves[$j]->{frequency} > $winAves[$j+1]->{frequency}) {
      } else {
        $failed = 1;
        last;
      }
    }
    next if (0!=$failed);
    for(my $j=$i-$windowFlank; $j<$i; ++$j) {
      if ($winAves[$j]->{frequency} < $winAves[$j+1]->{frequency}) {
      } else {
        $failed = 1;
        last;
      }
    }
    next if (0!=$failed);
    return $winAves[$i]->{kmercount};
  }
  return 0;
}

sub estimateErrorRate {
  my ($dataPointsRef) = @_;

  my $cov2CN = get2CNCoverage($dataPointsRef);
  my $cov1CN = ceil($cov2CN/2);
  my $covError = ceil($cov1CN/2);

  my $errorFrequency = 0;
  for(my $i=0; $i<$covError; ++$i) {
    if ($dataPointsRef->[$i]->{kmercount}<=$covError) {
      $errorFrequency += ($dataPointsRef->[$i]->{kmercount} * $dataPointsRef->[$i]->{frequency});
    } else {
      last;
    }
  }

  my $totalFrequency = 0;
  my $numDatapoints = scalar(@{$dataPointsRef});
  for(my $i=0; $i<$numDatapoints; ++$i) {
    $totalFrequency += ($dataPointsRef->[$i]->{kmercount} * $dataPointsRef->[$i]->{frequency});
  }

  if (0==$totalFrequency) {
    return -1;
  } else {
    return $errorFrequency / $totalFrequency;
  }
}
