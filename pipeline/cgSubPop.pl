#!/bin/env perl

##
## ecLego3 v3.2408.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##
## ecLego component
##

use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use FindBin;
use Storable qw(dclone); 

my $G_USAGE = "
$0 <command> -h
<command> [command-specific options]
setup    : set up scripts for subpopulation processing

      Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
      Please contact author for license (cheehongsg at gmail.com)
";

my $G_ECLEGO_PIPELINE_PATH = $FindBin::Bin;
my $G_ECLEGO_ROOT = dirname($G_ECLEGO_PIPELINE_PATH);
my $G_LOCAL_DIR = $G_ECLEGO_PIPELINE_PATH;

my $command = undef;
if (!defined $ARGV[0] || substr($ARGV[0],0,1) eq '-' || substr($ARGV[0],0,2) eq '--') {
	die("Please specify the command.\n",$G_USAGE);
}
$command = shift @ARGV; $command = lc $command;

# auto-flush for both stderr and stdout
select(STDERR);
$| = 1;
select(STDOUT);
$| = 1;

if ('setup' eq $command) {
	setupSubPopulations();
} else {
}

exit 0;

sub getAssemblies {
  my ($assembly, $subgraphIds, $rsgIds, $assembliesRef, $assembliesInfoRef) = @_;

  my %rsgIds = ();
  my @rsgIds = split(/,/,$rsgIds);
  for(my $i=0; $i<scalar(@rsgIds); ++$i) {
    if ($rsgIds[$i]=~/\=/) {
      my @bits = split(/\=/, $rsgIds[$i]);
      my @parts = split(/\@/, $bits[0]);
      my ($cgId,$rsgId) = ('','');
      # die "Invalid format for rsgId#",$i,"=",$rgsIds[$i],"\n"if (0==scalar(@parts));
      if (1==scalar(@parts)) {
        $rsgId = $parts[0];
        $cgId = $parts[0];
      } else {
        $rsgId = $parts[1];
        $cgId = $parts[0];
      }
      my @contigs = split(/;/,$bits[1]);
      grep { 
        if (!/^([+\-])*contig\_/) {
          my $firstChar = substr($_, 0 ,1);
          my $remaining = substr($_, 1);
          if ('+' eq $firstChar) {
            $_ = 'contig_'.$remaining;
          } elsif ('-' eq $firstChar) {
            $_ = '-contig_'.$remaining;
          } elsif ($firstChar =~ /[0..9]/) {
            $_ = 'contig_'.$_;
          } else {
            # pass thru as-is
          }
        }
      } @contigs;
      if (!exists $rsgIds{$rsgId}) {
        $rsgIds{$rsgId} = {order=>$i+1, contigs=>\@contigs, cgs=>[{id=>$cgId, order=>$i+1, contigs=>\@contigs}]};
      } else {
        # FIXME: ordering!
        push @{$rsgIds{$rsgId}->{cgs}}, {id=>$cgId, order=>$i+1, contigs=>\@contigs};
      }
    } else {
      $rsgIds{$rsgIds[$i]} = {order=>$i+1, cgs=>[]};
    }
  }
  if ($assembly=~/\.xls$/) {
    open INFILE, $assembly || die "Fail to open $assembly\n$!\n";
    my $header = 'subgraphId      numLinearVertices       numLinearEdges  numVertices     numEdges        num_0_chrom     kbp_0_chrom     avecov_0_chrom  num_1_low_amp   kbp_1_low_amp   avecov_1_low_amp      num_2_med_amp   kbp_2_med_amp   avecov_2_med_amp        num_3_high_amp  kbp_3_high_amp  avecov_3_high_amp       vertices        edges   cov_0_chrom     len_0_chrom   cov_1_low_amp   len_1_low_amp   cov_2_med_amp   len_2_med_amp   cov_3_high_amp  len_3_high_amp  linearVertices  linearEdges';
    my @headers = split(/\s+/, $header);
    my $numHeaders = scalar(@headers);
    while (<INFILE>) {
      next if (/^\/\//);
      chomp();
      my @bits = split(/\t/);
      my %item = ();
      for(my $i=0; $i<$numHeaders; ++$i) {
        $item{$headers[$i]} = $bits[$i];
      }
      if (exists $rsgIds{$item{subgraphId}}) {
        $item{sgId} = $subgraphIds;
        $item{rsgId} = $rsgIds{$item{subgraphId}}->{order};
        if (exists $rsgIds{$item{subgraphId}}->{contigs}) {
          # $item{contigs} = $rsgIds{$item{subgraphId}}->{contigs};
          # can be multiple!
          foreach my $cgsRef (@{$rsgIds{$item{subgraphId}}->{cgs}}) {
            my $itemRef = dclone(\%item);
            $itemRef->{contigs} = $cgsRef->{contigs};
            $itemRef->{assignedId} = $cgsRef->{id};
            push @{$assembliesRef}, $itemRef;
          }
        } else {
          my @contigs = split(/;/,$item{edges});
          grep { $_ = "contig_".$_ } @contigs;
          $item{contigs} = \@contigs;
          push @{$assembliesRef}, \%item;
        }
      }
    }
    close INFILE;
    @{$assembliesRef} = sort { $a->{sgId}<=>$b->{sgId} || $a->{rsgId}<=>$b->{rsgId} } @{$assembliesRef};

    my $tmpSeqFile = sprintf("%s/assembly.fasta", $assembliesInfoRef->{dirname});
    if (-f $tmpSeqFile) {
      $assembliesInfoRef->{file} = $tmpSeqFile;
    } else {
      my @bits = split(/\//, $assembliesInfoRef->{dirname});
      my $tmp = pop @bits;
      my $tmpDir = join("/", @bits);
      $tmpSeqFile = sprintf("%s/assembly.fasta", $tmpDir);
      if (-f $tmpSeqFile) {
        $assembliesInfoRef->{file} = $tmpSeqFile;
      } else {
        $assembliesInfoRef->{file} = '../assembly.fasta';
      }
    }
  } elsif ($assembly=~/\.fa$/ || $assembly=~/\.fasta$/) {
    open INFILE, $assembly || die "Fail to open $assembly\n$!\n";
    while (<INFILE>) {
      chomp();
      if (/^>/) {
        my @c=split(/\s+/);
        $c[0] =~ s/^>//;
        my %item=('contigs'=>[$c[0]]);
        while (my ($sgIds, $sgRef) = each %rsgIds) {
          my $found = 0;
          foreach my $ctg (@{$sgRef->{contigs}}) {
            if ($c[0] eq $ctg) {
              $item{sgId} = $sgIds;
              $item{rsgId} = $sgRef->{order};

              foreach my $cgRef (@{$sgRef->{cgs}}) {
                $found = 0;
                foreach my $cgCtg (@{$cgRef->{contigs}}) {
                  if ($c[0] eq $cgCtg) {
                    $item{assignedId} = $cgRef->{id};
                    $found = 1;
                    last;
                  }
                }
                last if (0!=$found);
              }

              $found = 1;
              last;
            }
          }
          last if (0!=$found);
        }
        push @{$assembliesRef}, \%item;
      }
    }
    close INFILE;

    $assembliesInfoRef->{file} = $assembly;
  } else {
    die "Unknown format $assembly\n";
  }
}

sub getExactContig {
  my ($contig) = @_;
  my $exactContig = $contig;
  if ($exactContig =~ /^[+\-]/) {
    $exactContig = substr($exactContig,1)
  }
  if ($exactContig =~ /:/) {
    my @bits = split(/:/, $exactContig);
    $exactContig = $bits[0];
  }
  return $exactContig;
}

sub retainAssembliesVariant {
  my ($vcf, $assembliesRef, $minsupport, $variantsRef) = @_;
  open INFILE, $vcf || die "Fail to open $vcf\n$!\n";

  $variantsRef->{headers} = ();
  $variantsRef->{rows} = ();
  my @cols = ();
  my $numCols = 0;
  my %contigs = ();
  foreach my $assemblyRef (@{$assembliesRef}) {
    foreach my $contig (@{$assemblyRef->{contigs}}) {
      my $exactContig = getExactContig($contig);
      if (!exists $contigs{$exactContig}) {
        $contigs{$exactContig} = [$assemblyRef];
      } else {
        # FIXME: if contig appears twice in the same assembly?
        push @{$contigs{$exactContig}}, $assemblyRef;
      }
    }
  }
  while (<INFILE>) {
    if (/^##/) {
      push @{$variantsRef->{headers}}, $_;
    } elsif (/^#/) {
      push @{$variantsRef->{headers}}, $_;
      my $header = $_;
      chomp($header);
      $header =~ s/^#//;
      @cols = split(/\t/, $header);
      $numCols = scalar(@cols);
      $variantsRef->{cols} = \@cols;
      $variantsRef->{numCols} = $numCols;
    } else {
      chomp();
      my @bits = split(/\t/);
      my %item = ();
      for(my $i=0; $i<$numCols; ++$i) {
        $item{$cols[$i]} = $bits[$i];
      }
      next if (! exists $contigs{$item{CHROM}});
      my ($numSupport) = $item{INFO} =~ /;SUPPORT=(\d+);/;
      next if ($numSupport<$minsupport);
      push @{$variantsRef->{rows}}, {line=>$_, itemRef=>\%item, assemblyRef=>$contigs{$item{CHROM}}, numSupport=>$numSupport};
    }
  }
  close INFILE;
}

sub writeCOIsVCF {
  my ($vcf, $assembliesRef, $minsupport, $variantsRef, $assembliesInfoRef) = @_;
  my %svtypesRows = ();
  foreach my $rowRef (@{$variantsRef->{rows}}) {
    my $variantRef = $rowRef->{itemRef};
    my ($svType) = $variantRef->{INFO} =~ /;SVTYPE=(\w+);/;
    if (!exists $svtypesRows{$svType}) {
      $svtypesRows{$svType} = [];
    }
    push @{$svtypesRows{$svType}}, $rowRef;
  }
  my @svtypesCount = ();
  while (my ($svtype, $rowsRef) = each %svtypesRows) {
    push @svtypesCount, {svtype=>$svtype, count=>scalar(@{$rowsRef})};
  }
  @svtypesCount = sort { $b->{count} <=> $a->{count}} @svtypesCount;

  my $vcfbasename = basename($vcf);
  if ($vcfbasename=~/\.vcf$/) {
    $vcfbasename =~ s/\.vcf$//;
  }
  # write the complete list
  $variantsRef->{summary} = [];
  my $ovcf = sprintf("%s/%s.ge%d.vcf", $assembliesInfoRef->{dirname}, $vcfbasename, $minsupport);
  my $basename = basename($ovcf);
  my $dirname = dirname($ovcf);
  my $line = sprintf "# VCF files\n# @ %s\n", $dirname;
  push @{$variantsRef->{summary}}, $line;
  $line = sprintf "# ALL support>=%d %d records %s\n", $minsupport, scalar(@{$variantsRef->{rows}}), $basename;
  push @{$variantsRef->{summary}}, $line;
  $variantsRef->{vcf} = $ovcf;
  open OUTFILE, ">$ovcf" || die "Fail to open $ovcf\n$!\n";
  grep { print OUTFILE $_; } @{$variantsRef->{headers}};
  grep { print OUTFILE $_->{line}, "\n"; } @{$variantsRef->{rows}};
  close OUTFILE;

  # write the svtype!
  foreach my $osvtypeRef (@svtypesCount) {
    $ovcf = sprintf("%s/%s.ge%d.%s.vcf", $assembliesInfoRef->{dirname}, $vcfbasename, $minsupport, $osvtypeRef->{svtype});
    $basename = basename($ovcf);
    $dirname = dirname($ovcf);
    $line = sprintf "# %s support>=%d %d records %s\n", $osvtypeRef->{svtype}, $minsupport, scalar(@{$svtypesRows{$osvtypeRef->{svtype}}}), $basename;
    push @{$variantsRef->{summary}}, $line;
    open OUTFILE, ">$ovcf" || die "Fail to open $ovcf\n$!\n";
    grep { print OUTFILE $_; } @{$variantsRef->{headers}};
    grep { print OUTFILE $_->{line}, "\n"; } @{$svtypesRows{$osvtypeRef->{svtype}}};
    close OUTFILE;
  }
  # print "\n";

  $variantsRef->{groups} = \%svtypesRows;
}

# TODO: test assemblies that require more than 1 contig!
sub writeCircularGenomes {
  my ($assembly, $assembliesRef, $variantsRef, $dicncov, $minsupport, $sample, $rotate, $svTypesRef, $minsvlen, $assembliesInfoRef) = @_;

  # $variantsRef->{groups} = \%svtypesRows;
  my $ofile = (0==$rotate) ? "setup.sh" : "setup.distilled.sh";
  $ofile = sprintf("%s/%s", $assembliesInfoRef->{dirname}, $ofile);
  open OUTFILE, ">$ofile" || die "Fail to open $ofile\n$!\n";

  # TODO: report header + statistics
  my $rotateTrial = 100000; # 100kb

  print OUTFILE "#####\n";
  print OUTFILE "# Circular Genomes generation\n";
  print OUTFILE "#####\n";

  print OUTFILE "\n";
  if (0==scalar(@{$variantsRef->{summary}})) {
    print OUTFILE "# There is no variant data available.\n";
  } else {
    print OUTFILE join("", @{$variantsRef->{summary}}), "\n";
  }


  print OUTFILE "\n";
  printf OUTFILE "export ECLEGO_PIPELINE=%s\n", $G_LOCAL_DIR;
  print OUTFILE "\n";

  print OUTFILE "\n";
  printf OUTFILE "export SIEVED_ASSEMBLY_FASTA=%s\n", $assembliesInfoRef->{file};
  my $seqFilePath = "\${SIEVED_ASSEMBLY_FASTA}";
  printf OUTFILE "export SIEVED_ASSEMBLY_PATH=%s\n", $assembliesInfoRef->{dirname};
  my $sievedAssemblyPath = "\${SIEVED_ASSEMBLY_PATH}";
  print OUTFILE "\n";

  for(my $i=0; $i<scalar(@{$assembliesRef}); ++$i) {
    my $assemblyRef = $assembliesRef->[$i];
    my $effectiveRSGID = (exists $assemblyRef->{assignedId}) ? $assemblyRef->{assignedId} : $assemblyRef->{rsgId};
    my $oprefix = sprintf("%s.sG%s_%s", $sample, $assemblyRef->{sgId}, $effectiveRSGID);
    my $cgCount = 0;

    my $fastaCollection = sprintf("%s.CGs.fa", $oprefix);
    my @collection = ();

    $cgCount++;
    my $cgName = '';
    my $activeContigs = scalar(@{$assemblyRef->{contigs}});
    my %activeContigs = ();
    # if ($activeContigs>1) {
    if ($activeContigs>0) {
      foreach my $contig (@{$assemblyRef->{contigs}}) {
        my $exactContig = getExactContig($contig);
        $cgName = sprintf("%s_%d_%s", , $oprefix, $cgCount, $exactContig);
        printf OUTFILE "# cG%d component: %s\n", $cgCount, $exactContig;
        print OUTFILE 'echo $(date) Generate CG ',$cgName, '..', "\n";
        print OUTFILE "perl \${ECLEGO_PIPELINE}/makeGenome.pl \\\n";
        printf OUTFILE "--seq %s \\\n", $seqFilePath;
        printf OUTFILE "--path \"%s\" \\\n", $exactContig;
        printf OUTFILE "--oprefix %s/%s\n", $sievedAssemblyPath, $cgName;
        print OUTFILE "\n";
        push @collection, sprintf("%s.fa", $cgName);
        push @collection, sprintf("%s.rc.fa", $cgName);

        my $activateContig = $contig;
        $activateContig =~ s/^\-//;
        $activeContigs{$activateContig}++;
      }
    }

    $cgName = sprintf("%s_%d", , $oprefix, $cgCount);
    printf OUTFILE "# make cG%d: %s\n", $cgCount, join(",", @{$assemblyRef->{contigs}});
    print OUTFILE 'echo $(date) Generate CG ',$cgName, '..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}/makeGenome.pl \\\n";
    printf OUTFILE "--seq %s \\\n", $seqFilePath;
    printf OUTFILE "--path \"%s\" \\\n", join(";", @{$assemblyRef->{contigs}});
    printf OUTFILE "--rotate %d \\\n", $rotate if (0!=$rotate);
    printf OUTFILE "--oprefix %s/%s\n", $sievedAssemblyPath, $cgName;
    print OUTFILE "\n";
    push @collection, sprintf("%s.fa", $cgName);
    push @collection, sprintf("%s.rc.fa", $cgName);

    # rotate backbone by 100kb
    $cgName = sprintf("%s_%d_rot100k", , $oprefix, $cgCount);
    printf OUTFILE "# make cG%d: %s\n", $cgCount, join(",", @{$assemblyRef->{contigs}});
    print OUTFILE 'echo $(date) Generate CG ',$cgName, '..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}/makeGenome.pl \\\n";
    printf OUTFILE "--seq %s \\\n", $seqFilePath;
    printf OUTFILE "--path \"%s\" \\\n", join(";", @{$assemblyRef->{contigs}});
    printf OUTFILE "--rotate %d \\\n", $rotateTrial + $rotate;
    printf OUTFILE "--oprefix %s/%s\n", $sievedAssemblyPath, $cgName;
    print OUTFILE "\n";
    push @collection, sprintf("%s.fa", $cgName);
    push @collection, sprintf("%s.rc.fa", $cgName);

    # enumerate the SVs
    # my @svTypes = ('INV', 'DEL', 'DUP');
    # my @svTypes = ('INV', 'DEL');
    my @svTypes = (); push @svTypes, @{$svTypesRef};
    foreach my $svType (@svTypes) {
      next if (!exists $variantsRef->{groups}->{$svType});

      # decreasing copy numbers
      my $svVariantsRef = $variantsRef->{groups}->{$svType};
      my @orders = sort { $b->{numSupport}<=>$a->{numSupport} } @{$svVariantsRef};
      foreach my $ovarRef (@orders) {
        my $itemRef = $ovarRef->{itemRef};
        if ($activeContigs>0) {
          next if (! exists $activeContigs{$itemRef->{CHROM}});
        }

        my ($svlen) = $itemRef->{INFO} =~ /;SVLEN=([\-+]*\d+);/;
        next if (($svlen<0?-1*$svlen:$svlen)<$minsvlen);
        
        my ($endpos) = $itemRef->{INFO} =~ /;END=(\d+);/;
        $cgCount++;
        $cgName = sprintf("%s_%d", , $oprefix, $cgCount);
        printf OUTFILE "# make cG%d: %s + %s = %s:%d-%d SVLEN=%d SUPPORT=%d (%.1fc)\n", 
          $cgCount, join(",", @{$assemblyRef->{contigs}}),
          $itemRef->{ID}, $itemRef->{CHROM}, $itemRef->{POS}, $endpos, $svlen, $ovarRef->{numSupport}, $ovarRef->{numSupport}/$dicncov*2;
        my $cgDetails = sprintf("%s + %s = %s:%d-%d SVLEN=%d SUPPORT=%d (%.1fc)", join(",", @{$assemblyRef->{contigs}}),
          $itemRef->{ID}, $itemRef->{CHROM}, $itemRef->{POS}, $endpos, $svlen, $ovarRef->{numSupport}, $ovarRef->{numSupport}/$dicncov*2);
        print OUTFILE 'echo $(date) Generate CG ',$cgName, '.. \'', $cgDetails, '\'', "\n";
        # FIXME: should figure out the impact of large deletion w.r.t. rotation!
        my $toRotate = (0!=$rotate) ? 1 : 0;
        my $adjustment = 0;
        if ($itemRef->{ID}=~/\.DEL\./i) {
          if (abs($svlen)>$rotate) {
            $toRotate = 0;
            printf OUTFILE "echo WARNING: svlen=%d \\> rotation=%d, skipping rotation.\n", abs($svlen), $rotate;
          } else {
            $adjustment = $svlen;
          }
        } elsif ($itemRef->{ID}=~/\.INV\./i) {

        } else {
          # don't rotate
          $toRotate = 0;
        }
        print OUTFILE "perl \${ECLEGO_PIPELINE}/makeGenome.pl \\\n";
        printf OUTFILE "--seq %s \\\n", $seqFilePath;
        my $vcfbasename = basename($variantsRef->{vcf});
        my $vcfdirname = dirname($variantsRef->{vcf});
        # printf OUTFILE "--vcf %s \\\n", $variantsRef->{vcf};
        printf OUTFILE "--vcf %s/%s \\\n", $sievedAssemblyPath, $vcfbasename;
        my @paths = ();
        foreach my $contig (@{$assemblyRef->{contigs}}) {
          my $norm_contig = $contig; $norm_contig =~ s/^\-//;
          if ($itemRef->{CHROM} eq $norm_contig) {
            push @paths, sprintf("%s:%s", $contig, $itemRef->{ID});
          } else {
            push @paths, $contig;
          }
        }
        printf OUTFILE "--path \"%s\" \\\n", join(";", @paths);
        printf OUTFILE "--rotate %d \\\n", ($rotate + $adjustment) if (0!=$toRotate);
        printf OUTFILE "--oprefix %s/%s\n", $sievedAssemblyPath, $cgName;
        print OUTFILE "\n";
        push @collection, sprintf("%s.fa", $cgName);
        push @collection, sprintf("%s.rc.fa", $cgName);
        print OUTFILE "\n";
      }
    }

    printf OUTFILE "# collate all CGs of %s (%d files)\n", $oprefix, scalar(@collection);
    print OUTFILE 'echo $(date) Generate ',$fastaCollection, '..', "\n";
    # print OUTFILE "cat ", join(" ", @collection), " > ", $fastaCollection, "\n";
    for(my $i=0; $i<scalar(@collection); ++$i) {
      if (0==$i) {
        print OUTFILE "cat $sievedAssemblyPath/$collection[$i] > $sievedAssemblyPath/$fastaCollection\n";
      } else {
        print OUTFILE "cat $sievedAssemblyPath/$collection[$i] >> $sievedAssemblyPath/$fastaCollection\n";
      }
    } 
    print OUTFILE "\n";
  }

  print OUTFILE 'echo $(date) Done.', "\n";
  print OUTFILE "\n";

  close OUTFILE;

  # if (0==$rotate) {
  #   print "./setup.sh 2>&1 | tee setup.run.log\n"
  # } else {
  #   print "./setup.distilled.sh 2>&1 | tee setup.distilled.run.log\n"
  # }
}


sub writeMapSVCall {
  my ($assembliesRef, $sample, $rotate, $assembliesInfoRef, $genome) = @_;

  my $ofile = (0==$rotate) ? "process.sh" : "process.sh";
  $ofile = sprintf("%s/%s", $assembliesInfoRef->{dirname}, $ofile);
  open OUTFILE, ">$ofile" || die "Fail to open $ofile\n$!\n";
  # print "\n";
  # print "\n";

  print OUTFILE "#####\n";
  print OUTFILE "# Circular Genomes mapping and structural variations calling\n";
  print OUTFILE "#####\n";

  print OUTFILE "\n";
  printf OUTFILE "export ECLEGO_PIPELINE_PATH=%s\n", $G_LOCAL_DIR;
  print OUTFILE 'export ECLEGO_ROOT=$(dirname "${ECLEGO_PIPELINE_PATH}")', "\n";
  printf OUTFILE "export REFGENOME=%s\n", $genome;
  print OUTFILE "\n";
  print OUTFILE 'if command -v apptainer &> /dev/null; then', "\n";
  print OUTFILE '    export APPTAINER_BINDPATH=${ECLEGO_ROOT}', "\n";
  print OUTFILE 'else', "\n";
  print OUTFILE '    export SINGULARITY_BINDPATH=${ECLEGO_ROOT}', "\n";
  print OUTFILE 'fi', "\n";
  print OUTFILE "\n";
  print OUTFILE 'export MINIMAP2CMD="${ECLEGO_PIPELINE_PATH}/minimap2_2.24.sif minimap2"', "\n";
  print OUTFILE 'export MINIMAP2THREADS=4', "\n";
  print OUTFILE 'export SAMTOOLSCMD="${ECLEGO_PIPELINE_PATH}/samtools_1.19.sif samtools"', "\n";
  print OUTFILE 'export SAMTOOLSTHREADS=4', "\n";
  print OUTFILE 'export SNIFFLESSCMD="${ECLEGO_PIPELINE_PATH}/sniffles_2.3.3.sif sniffles"', "\n";

  print OUTFILE "\n";
  printf OUTFILE "export SIEVED_ASSEMBLY_FASTA=%s\n", $assembliesInfoRef->{file};
  my $seqFilePath = "\${SIEVED_ASSEMBLY_FASTA}";
  printf OUTFILE "export SIEVED_ASSEMBLY_PATH=%s\n", $assembliesInfoRef->{dirname};
  my $sievedAssemblyPath = "\${SIEVED_ASSEMBLY_PATH}";
  print OUTFILE "\n";

  for(my $i=0; $i<scalar(@{$assembliesRef}); ++$i) {
    my $assemblyRef = $assembliesRef->[$i];
    my $effectiveRSGID = (exists $assemblyRef->{assignedId}) ? $assemblyRef->{assignedId} : $assemblyRef->{rsgId};
    my $oprefix = sprintf("%s.sG%s_%s", $sample, $assemblyRef->{sgId}, $effectiveRSGID);
    my $fastaCollection = sprintf("%s.CGs.fa", $oprefix);
    my $bamFile = sprintf("%s.CGs.bam", $oprefix);

    print OUTFILE "\n";
    print OUTFILE "# ",$fastaCollection, "\n";
    printf OUTFILE "export FASTA=%s\n", $fastaCollection;
    print OUTFILE "\n";

    print OUTFILE "# mapping\n";
    print OUTFILE 'echo $(date) Mapping ',$fastaCollection, ' against T2Tv2..', "\n";
    print OUTFILE '${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} -L "${REFGENOME/%.fasta/.mmi}" ${SIEVED_ASSEMBLY_PATH}/${FASTA} \\', "\n";
    print OUTFILE '| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${SIEVED_ASSEMBLY_PATH}/${FASTA/%.fa/.bam} ; \\', "\n";
    print OUTFILE '${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${SIEVED_ASSEMBLY_PATH}/${FASTA/%.fa/.bam}', "\n";
    print OUTFILE "\n";

    print OUTFILE "# sv calling\n";
    printf OUTFILE "export BAMFILE=%s\n", $bamFile;
    print OUTFILE 'echo $(date) SV Calling pass minsupport=1..', "\n";
    print OUTFILE '${SNIFFLESSCMD} --allow-overwrite --minsupport 1 --reference ${REFGENOME} --input ${SIEVED_ASSEMBLY_PATH}/${BAMFILE} --vcf ${SIEVED_ASSEMBLY_PATH}/${BAMFILE}.pass.vcf', "\n";
    print OUTFILE "\n";

    print OUTFILE "# sv calling for various lengths\n";
    print OUTFILE 'for i in 20 10 5 1 ; do \\', "\n";
    print OUTFILE 'echo $(date) SV Calling pass minsupport=1 minsvlen=${i}.. ; \\', "\n";
    print OUTFILE '${SNIFFLESSCMD} --allow-overwrite --minsupport 1 --minsvlen ${i} --reference ${REFGENOME} --input ${SIEVED_ASSEMBLY_PATH}/${BAMFILE} --vcf ${SIEVED_ASSEMBLY_PATH}/${BAMFILE}.pass.minlen${i}.vcf ; \\', "\n";
    print OUTFILE 'done', "\n";     
    print OUTFILE "\n";

    print OUTFILE "# for comprehensiveness\n";
    print OUTFILE '${SNIFFLESSCMD} --qc-output-all --allow-overwrite --reference ${REFGENOME} --input ${SIEVED_ASSEMBLY_PATH}/${BAMFILE} --vcf ${SIEVED_ASSEMBLY_PATH}/${BAMFILE}.raw.vcf', "\n";
    print OUTFILE "\n";

    print OUTFILE "# for support read ids (future use))\n";
    print OUTFILE '${SNIFFLESSCMD} --output-rnames --qc-output-all --allow-overwrite --minsupport 1 --reference ${REFGENOME} --input ${SIEVED_ASSEMBLY_PATH}/${BAMFILE} --vcf ${SIEVED_ASSEMBLY_PATH}/${BAMFILE}.pass.vcf.temp', "\n";
    print OUTFILE "\n";
  }

  print OUTFILE 'echo $(date) Done.', "\n";
  print OUTFILE "\n";

  close OUTFILE;
}


sub writePublishForReview {
  my ($assembliesRef, $sample, $rotate, $assembliesInfoRef) = @_;

  my $ofile = (0==$rotate) ? "process.sh" : "process.sh";
  $ofile = sprintf("%s/%s", $assembliesInfoRef->{dirname}, $ofile);
  open OUTFILE, ">>$ofile" || die "Fail to open $ofile\n$!\n";
  print OUTFILE "\n";
  print OUTFILE "\n";
  print OUTFILE "#####\n";
  print OUTFILE "# Circular Genomes baseline for review\n";
  print OUTFILE "#####\n";

  print OUTFILE "\n";
  printf OUTFILE "export ECLEGO_PIPELINE=%s\n", $G_LOCAL_DIR;
  
  print OUTFILE "\n";
  printf OUTFILE "export SIEVED_ASSEMBLY_FASTA=%s\n", $assembliesInfoRef->{file};
  my $seqFilePath = "\${SIEVED_ASSEMBLY_FASTA}";
  printf OUTFILE "export SIEVED_ASSEMBLY_PATH=%s\n", $assembliesInfoRef->{dirname};
  my $sievedAssemblyPath = "\${SIEVED_ASSEMBLY_PATH}";
  print OUTFILE "\n";

  for(my $i=0; $i<scalar(@{$assembliesRef}); ++$i) {
    my $assemblyRef = $assembliesRef->[$i];
    my $effectiveRSGID = (exists $assemblyRef->{assignedId}) ? $assemblyRef->{assignedId} : $assemblyRef->{rsgId};
    my $cgPrefix = sprintf("%s.sG%s_%s.CGs", $sample, $assemblyRef->{sgId}, $effectiveRSGID);
    my $fastaCollection = sprintf("%s.fa", $cgPrefix);
    my $bamFile = sprintf("%s.bam", $cgPrefix);

    print OUTFILE "#\n";
    print OUTFILE "# ",$bamFile, " from ", $fastaCollection, "\n";
    print OUTFILE "#\n";
    print OUTFILE "\n";

    print OUTFILE "# all alignments\n";
    print OUTFILE 'echo $(date) Generate All alignments table..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}cgbam.pl tablebam \\\n";
    printf OUTFILE "--bam %s/%s \\\n",
      $sievedAssemblyPath, $bamFile;
    printf OUTFILE "> %s/%s.alignment.table.all.xls\n",
      $sievedAssemblyPath, $cgPrefix;
    print OUTFILE "\n";

    print OUTFILE "# breakpoints\n";
    print OUTFILE 'echo $(date) Generate breakpoints..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}cgbam.pl breakpoints \\\n";
    printf OUTFILE "--alignment %s/%s.alignment.table.all.xls --bed \\\n",
      $sievedAssemblyPath, $cgPrefix;
    printf OUTFILE "| tee %s/%s.breakpoints.xls\n",
      $sievedAssemblyPath, $cgPrefix;
    print OUTFILE "\n";

    print OUTFILE 'echo $(date) Done.', "\n";
    print OUTFILE "\n";
  }

  close OUTFILE;



  ##### iterative refinement process

  $ofile = (0==$rotate) ? "refine.sh" : "refine.distilled.sh";
  $ofile = sprintf("%s/%s", $assembliesInfoRef->{dirname}, $ofile);
  open OUTFILE, ">$ofile" || die "Fail to open $ofile\n$!\n";
  # print "\n";
  # print "\n";
  print OUTFILE "#####\n";
  print OUTFILE "# Circular Genomes refinements\n";
  print OUTFILE "#####\n";
    
  print OUTFILE "\n";
  printf OUTFILE "export ECLEGO_PIPELINE=%s\n", $G_LOCAL_DIR;
  
  print OUTFILE "\n";
  printf OUTFILE "export SIEVED_ASSEMBLY_FASTA=%s\n", $assembliesInfoRef->{file};
  printf OUTFILE "export SIEVED_ASSEMBLY_PATH=%s\n", $assembliesInfoRef->{dirname};
  print OUTFILE "\n";

  for(my $i=0; $i<scalar(@{$assembliesRef}); ++$i) {
    my $assemblyRef = $assembliesRef->[$i];
    my $effectiveRSGID = (exists $assemblyRef->{assignedId}) ? $assemblyRef->{assignedId} : $assemblyRef->{rsgId};
    my $cgPrefix = sprintf("%s.sG%s_%s.CGs", $sample, $assemblyRef->{sgId}, $effectiveRSGID);
    my $fastaCollection = sprintf("%s.fa", $cgPrefix);
    my $bamFile = sprintf("%s.bam", $cgPrefix);

    print OUTFILE "#\n";
    print OUTFILE "# ",$bamFile, " from ", $fastaCollection, "\n";
    print OUTFILE "#\n";
    print OUTFILE "\n";

    print OUTFILE "# preferred (start with + strand) alignments\n";
    printf OUTFILE "# --exclude=%s.sG%s_%s_1_rot100k --exclude=%s.sG%s_%s_1_rot100k/rc \\\n", 
      $sample, $assemblyRef->{sgId}, $effectiveRSGID,
      $sample, $assemblyRef->{sgId}, $effectiveRSGID;
    print OUTFILE 'echo $(date) Generate strand selected alignments table..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}cgbam.pl tablebam \\\n";
    printf OUTFILE "--bam %s/%s --easycoord \\\n", $sievedAssemblyPath, $bamFile;
    printf OUTFILE "--exclude=%s.sG%s_%s_1_rot100k --exclude=%s.sG%s_%s_1_rot100k/rc \\\n", 
      $sample, $assemblyRef->{sgId}, $effectiveRSGID,
      $sample, $assemblyRef->{sgId}, $effectiveRSGID
      if (0!=$rotate);
    printf OUTFILE "> %s/%s.alignment.table.xls\n", $sievedAssemblyPath, $cgPrefix;
    print OUTFILE "\n";
    
    my $minSVLen = 3000;
    my $minSVLenLabel = '3K';
    printf OUTFILE "# breakpoints with svlen>=%d-bp\n", $minSVLen;
    print OUTFILE 'echo $(date) Generate breakpoints with svlen\>=',$minSVLen,'-bp..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}cgbam.pl breakpoints \\\n";
    printf OUTFILE "--alignment %s/%s.alignment.table.xls \\\n", $sievedAssemblyPath, $cgPrefix;
    printf OUTFILE "--bed --mindellen %d --mininslen %d \\\n", $minSVLen, $minSVLen;
    printf OUTFILE "| tee %s/%s.breakpoints.delge%s.insge%s.xls\n",
      $sievedAssemblyPath, $cgPrefix, $minSVLenLabel, $minSVLenLabel;
    printf OUTFILE "mv %s/%s.alignment.table.xls.bed %s/%s.breakpoints.delge%s.insge%s.bed\n",
      $sievedAssemblyPath, $cgPrefix, $sievedAssemblyPath, $cgPrefix, $minSVLenLabel, $minSVLenLabel;
    print OUTFILE "\n";

    $minSVLen = 300;
    $minSVLenLabel = '300';
    printf OUTFILE "# breakpoints with svlen>=%d-bp\n", $minSVLen;
    print OUTFILE 'echo $(date) Generate breakpoints with svlen\>=',$minSVLen,'-bp..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}cgbam.pl breakpoints \\\n";
    printf OUTFILE "--alignment %s/%s.alignment.table.xls \\\n", $sievedAssemblyPath, $cgPrefix;
    printf OUTFILE "--bed --mindellen %d --mininslen %d \\\n", $minSVLen, $minSVLen;
    printf OUTFILE "| tee %s/%s.breakpoints.delge%s.insge%s.xls\n", 
      $sievedAssemblyPath, $cgPrefix, $minSVLenLabel, $minSVLenLabel;
    printf OUTFILE "mv %s/%s.alignment.table.xls.bed %s/%s.breakpoints.delge%s.insge%s.bed\n",
      $sievedAssemblyPath, $cgPrefix, $sievedAssemblyPath, $cgPrefix, $minSVLenLabel, $minSVLenLabel;
    print OUTFILE "\n";

    $minSVLen = 300;
    $minSVLenLabel = '300';
    printf OUTFILE "# using breakpoints with svlen>=%d-bp as drill down set\n", $minSVLen;
    print OUTFILE 'echo $(date) Generate final breakpoints with svlen\>=',$minSVLen,'-bp..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}cgbam.pl breakpoints \\\n";
    printf OUTFILE "--alignment %s/%s.alignment.table.xls \\\n", $sievedAssemblyPath, $cgPrefix;
    printf OUTFILE "--bed --mindellen %d --mininslen %d \\\n", $minSVLen, $minSVLen;
    printf OUTFILE "| tee %s/%s.breakpoints.xls\n", $sievedAssemblyPath, $cgPrefix;
    printf OUTFILE "mv %s/%s.alignment.table.xls.bed %s/%s.breakpoints.bed\n",
      $sievedAssemblyPath, $cgPrefix, $sievedAssemblyPath, $cgPrefix;
    print OUTFILE "\n";

    print OUTFILE "# tag the alignment with auxiliary info\n";
    print OUTFILE 'echo $(date) Tagging bam file for visualization..', "\n";
    print OUTFILE "perl \${ECLEGO_PIPELINE}cgbam.pl tagbam \\\n";
    printf OUTFILE "--bam %s/%s \\\n", $sievedAssemblyPath, $bamFile;
    printf OUTFILE "--alignment %s/%s.alignment.table.xls\n", $sievedAssemblyPath, $cgPrefix;
    print OUTFILE "\n";
    
    print OUTFILE 'echo $(date) Done.', "\n";
    print OUTFILE "\n";
  }

  close OUTFILE;
}


sub setupSubPopulations {
  my $G_USAGE = "
$0 download --assembly <bamFile> --easycoord --bed
  --assembly STR      the assembly results
  --vcf STR           the variants call file
  --subgraph STR      the subgraph id
  --rsg STR           the refined subgraph id
  --dicncov INT       read coverage per 2 copy number
  
  --minsupport INT    min number of supporting read
  --sample STR        sample id
  --scriptPath STR    path to script
  --rotate INT        generate distilled set with rotation
";

  my $assembly = undef;
  my $vcf = undef;
  my $subgraphIds = undef;
  my $rsgIds = undef;
  my $sample = undef;
  my $scriptPath = undef;
  my $dicncov = 0;
  my $genome = undef;
  my $minsupport = 0;
  my $rotate = 0;
  my $minsvlen = 0;
  my @svTypes = ();
  my $verbose = 0;
	my $help = 0;
	
	GetOptions (
	"assembly=s"   => \$assembly,
	"vcf=s"        => \$vcf,
	"subgraph=s"   => \$subgraphIds,
	"rsg=s"        => \$rsgIds,
	"sample=s"     => \$sample,
	"scriptPath=s" => \$scriptPath,
	"dicncov=i"    => \$dicncov,
	"genome=s"     => \$genome,
	"minsupport=i" => \$minsupport,
	"rotate=i"     => \$rotate,
	"sv=s"         => \@svTypes,
	"minsvlen=i"   => \$minsvlen,
	"verbose!"     => \$verbose,
	"help!"        => \$help)
	or die("Error in command line arguments\n$G_USAGE");
	
	die "$G_USAGE" if ($help);
  die "coverage must be a positive number: $dicncov\n" if ($dicncov<1);
  if (0==$minsupport) {
    $minsupport = int($dicncov/2); # minimumly one copy number
    if ($minsupport<4) {
      $minsupport = 4;
    }
  }
  if (defined $scriptPath) {
    die "Path does not exist! $scriptPath\n" if (! -d $scriptPath);
    if (substr($scriptPath,-1) ne '/') {
      $scriptPath .= '/';
    }
    $G_LOCAL_DIR = $scriptPath;
  }
  if (!defined $genome) {
    $genome = sprintf("%s%s", $G_LOCAL_DIR, 'genomes/t2tv2.fasta');
  }
  if (! -f $genome) {
    die "Genome file does not exist! $genome\n";
  }
  push @svTypes, 'INV', 'DEL' if (0==scalar(@svTypes));

  my $basename = basename($assembly);
  my $dirname = dirname($assembly);
  my $logfile = $basename;
  $logfile .= ".process.log";

  # generate the vcf files (svttype)
  # generate (scripts) the table of genomes generated for cross-check
  # separate script to map genomes back to t2tv2
  my @assemblies = ();
  my %assembliesInfo = (dirname=>$dirname);
  getAssemblies($assembly, $subgraphIds, $rsgIds, \@assemblies, \%assembliesInfo);
  my %variants = ();
  retainAssembliesVariant($vcf, \@assemblies, $minsupport, \%variants);

  # write COIs VCF
  # write COIs SVTYPE VCF
  writeCOIsVCF($vcf, \@assemblies, $minsupport, \%variants, \%assembliesInfo);

  # write CGs
  writeCircularGenomes($assembly, \@assemblies, \%variants, $dicncov, $minsupport, $sample, $rotate, \@svTypes, $minsvlen, \%assembliesInfo);

  # write downstream processing
  writeMapSVCall(\@assemblies, $sample, $rotate, \%assembliesInfo, $genome);

  # write downstream publishing
  writePublishForReview(\@assemblies, $sample, $rotate, \%assembliesInfo);
}
