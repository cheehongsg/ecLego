##
## ecLego2 v2.04.01
## Copyright (C) 2019-2023 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##
## ecLego component
##

use strict;
use Data::Dumper;
use POSIX qw/ceil/;

my $kmerGraphFile = $ARGV[0];
open INFILE, $kmerGraphFile || die "Fail to open $kmerGraphFile\n$!\n";
my @counts = ();
while (<INFILE>) {
  chomp();
  my ($kmerCount, $kmerFrequency) = split(/\t/);
  push @counts, {kmercount=>$kmerCount, frequency=>$kmerFrequency};
 
}
close INFILE;

my $cov2CN = get2CNCoverage(\@counts);
my $cov1CN = ceil($cov2CN/2);
my $cov4CN = ceil($cov2CN*2);
my $cov5CN = ceil($cov2CN*5/2);
my $cov6CN = ceil($cov2CN*3);
my $errorRate = estimateErrorRate(\@counts);

printf "# %s\n", $kmerGraphFile;
printf "covCN=%d\n", $cov1CN;
printf "cov2CN=%d\n", $cov2CN;
printf "cov4CN=%d\n", $cov4CN;
printf "cov5CN=%d\n", $cov5CN;
printf "cov6CN=%d\n", $cov6CN;
printf "errorRate=%.3f\n", $errorRate;

exit 0;

sub get2CNCoverage {
  my ($dataPointsRef) = @_;

  my $numDatapoints = scalar(@{$dataPointsRef});
  for(my $i=1; $i<($numDatapoints-1); ++$i) {
    if (($dataPointsRef->[$i-1]->{frequency} < $dataPointsRef->[$i]->{frequency})
      && ($dataPointsRef->[$i]->{frequency} > $dataPointsRef->[$i+1]->{frequency})) {
      return $dataPointsRef->[$i]->{kmercount};
    }
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
