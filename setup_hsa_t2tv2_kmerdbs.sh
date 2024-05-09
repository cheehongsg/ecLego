#!/bin/bash
set -eu -o pipefail

##
## ecLego2 v2.04.01
## Copyright (C) 2019-2023 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

echo $(date) - Downloading Human T2Tv2 annotation genomes from UCSC..
curl -o pipeline/resources/censat.bb "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.censat/censat.bb" -C -

echo $(date) - Renaming chr\{1..22\} to \{1..22\}, chr\{X,Y,M\} to \{X,Y,MT\}..
export SINGULARITY_BINDPATH=$(pwd)
./pipeline/ucsc-bigbedtobed_377.sif bigBedToBed pipeline/resources/censat.bb pipeline/resources/censat.bed.tmp
sed s?^chr?? pipeline/resources/censat.bed.tmp > pipeline/resources/censat.bed
rm pipeline/resources/censat.bed.tmp

echo $(date) - Generate greylist..
#
cat pipeline/resources/censat.bed \
| perl -ne 'BEGIN{$MINSPAN=20000;}{@c=split(/\t/);
    if ($c[3]=~/^hor/ || $c[3]=~/^hsat[23]\_/ || $c[3]=~/^hsat1[AB]\_/ || $c[3]=~/^rDNA\_/) { print $_; } 
    elsif ($c[3]=~/^censat\_/ || $c[3]=~/^dhor\_/ || $c[3]=~/^mon\_/) { print $_ if (($c[2]-$c[1])>=20000); }
}' \
| ./pipeline/bedtools_2.31.0.sif sortBed -i - \
| ./pipeline/bedtools_2.31.0.sif mergeBed -i - \
> pipeline/resources/t2tv2.greylist.bed
#
./pipeline/seqtk_1.3.sif seqtk subseq \
pipeline/genomes/t2tv2.fasta \
pipeline/resources/t2tv2.greylist.bed \
> pipeline/resources/t2tv2.greylist.fasta

echo $(date) - Generate greylist kmer database..
mkdir -p pipeline/resources/t2tv2.greylist.tmp/
./pipeline/kmc_3.2.1.sif kmc -k29 -m12 -t12 -ci1 -cs1000000 \
-fa pipeline/resources/t2tv2.greylist.fasta \
pipeline/resources/t2tv2.greylist.k29 pipeline/resources/t2tv2.greylist.tmp/

echo $(date) - Generate genome kmer database..
mkdir -p pipeline/resources/t2tv2.tmp/
./pipeline/kmc_3.2.1.sif kmc -k29 -m12 -t12 -ci1 -cs1000000 \
-fm pipeline/genomes/t2tv2.fasta \
pipeline/resources/t2tv2.k29 pipeline/resources/t2tv2.tmp/

echo $(date) - End of kmers databases setup.

