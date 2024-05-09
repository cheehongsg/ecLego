#!/bin/bash
set -eu -o pipefail

##
## ecLego2 v2.04.01
## Copyright (C) 2019-2023 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

echo $(date) - Downloading Human T2Tv2 genomes..
curl -o pipeline/genomes/chm13v2.0.fa.original.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -C -

echo $(date) - Renaming chr\{1..22\} to \{1..22\}, chr\{X,Y,M\} to \{X,Y,MT\}..
gzip -dc pipeline/genomes/chm13v2.0.fa.original.gz \
| perl -ne 'if (/^>/) { chomp(); ($_) = split(/\s+/); if (/>chrM/) { $_=~s/>chrM/>MT/; } else { $_=~s/>chr/>/; } $_.="\n"; } print $_;' \
> pipeline/genomes/t2tv2.fasta

export SINGULARITY_BINDPATH=$(pwd)
echo $(date) - Indexing Human T2Tv2 genomes..
./pipeline/samtools_v1.15.1.sif samtools faidx pipeline/genomes/t2tv2.fasta

echo $(date) - End of genome setup.
