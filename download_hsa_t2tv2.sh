#!/bin/bash
set -eu -o pipefail

##
## ecLego3 v3.2404.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

export ECGENOME_FOLDER=pipeline/genomes
[ ! -d "${ECGENOME_FOLDER}" ] && mkdir -p "${ECGENOME_FOLDER}"

echo $(date) - Downloading Human T2Tv2 genomes..
curl -# -L -o ${ECGENOME_FOLDER}/chm13v2.0.fa.original.gz -w "%{filename_effective}\n" https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -C -

echo $(date) - Renaming chr\{1..22\} to \{1..22\}, chr\{X,Y,M\} to \{X,Y,MT\}..
gzip -dc ${ECGENOME_FOLDER}/chm13v2.0.fa.original.gz \
| perl -ne 'if (/^>/) { chomp(); ($_) = split(/\s+/); if (/>chrM/) { $_=~s/>chrM/>MT/; } else { $_=~s/>chr/>/; } $_.="\n"; } print $_;' \
> ${ECGENOME_FOLDER}/t2tv2.fasta

export SINGULARITY_BINDPATH=$(pwd)
export APPTAINER_BINDPATH=${SINGULARITY_BINDPATH}
echo $(date) - Indexing Human T2Tv2 genomes..
./pipeline/samtools_1.19.sif samtools faidx ${ECGENOME_FOLDER}/t2tv2.fasta

echo $(date) - Downloading Human MT genomes..
curl -# -L -o ${ECGENOME_FOLDER}/MT.fasta -w "%{filename_effective}\n" "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_012920.1&rettype=fasta&retmode=text"

echo $(date) - Indexing Human MT genomes..
./pipeline/samtools_1.19.sif samtools faidx ${ECGENOME_FOLDER}/MT.fasta

echo $(date) - End of genome setup.
