#!/bin/bash
set -eu -o pipefail

##
## ecLego3 v3.2404.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

export ECGENOME_FOLDER=pipeline/genomes
export ECRESOURCE_FOLDER=pipeline/resources
[ ! -d "${ECRESOURCE_FOLDER}" ] && mkdir -p "${ECRESOURCE_FOLDER}"

echo $(date) - Downloading Human T2Tv2 annotation genomes from UCSC..
curl -# -L -o ${ECRESOURCE_FOLDER}/censat.bb -w "%{filename_effective}\n" "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.censat/censat.bb" -C -

echo $(date) - Renaming chr\{1..22\} to \{1..22\}, chr\{X,Y,M\} to \{X,Y,MT\}..
if command -v apptainer &> /dev/null; then
    export APPTAINER_BINDPATH=$(pwd)
else
    export SINGULARITY_BINDPATH=$(pwd)
fi
./pipeline/ucsc-bigbedtobed_377.sif bigBedToBed ${ECRESOURCE_FOLDER}/censat.bb ${ECRESOURCE_FOLDER}/censat.bed.tmp
sed s?^chr?? ${ECRESOURCE_FOLDER}/censat.bed.tmp > ${ECRESOURCE_FOLDER}/censat.bed
rm ${ECRESOURCE_FOLDER}/censat.bed.tmp

echo $(date) - Generate greylist..
#
cat ${ECRESOURCE_FOLDER}/censat.bed \
| perl -ne 'BEGIN{$MINSPAN=20000;}{@c=split(/\t/);
    if ($c[3]=~/^hor/ || $c[3]=~/^hsat[23]\_/ || $c[3]=~/^hsat1[AB]\_/ || $c[3]=~/^rDNA\_/) { print $_; } 
    elsif ($c[3]=~/^censat\_/ || $c[3]=~/^dhor\_/ || $c[3]=~/^mon\_/) { print $_ if (($c[2]-$c[1])>=20000); }
}' \
| ./pipeline/bedtools_2.31.0.sif sortBed -i - \
| ./pipeline/bedtools_2.31.0.sif mergeBed -i - \
> ${ECRESOURCE_FOLDER}/t2tv2.greylist.bed
#
./pipeline/seqtk_1.3.sif seqtk subseq \
${ECGENOME_FOLDER}/t2tv2.fasta \
${ECRESOURCE_FOLDER}/t2tv2.greylist.bed \
> ${ECRESOURCE_FOLDER}/t2tv2.greylist.fasta

echo $(date) - Generate greylist kmer database..
export ECTMP_FOLDER=${ECRESOURCE_FOLDER}/t2tv2.greylist.tmp/
[ ! -d "${ECTMP_FOLDER}" ] && mkdir -p "${ECTMP_FOLDER}"
./pipeline/kmc_3.2.1.sif kmc -k29 -m12 -t12 -ci1 -cs1000000 \
-fa ${ECRESOURCE_FOLDER}/t2tv2.greylist.fasta \
${ECRESOURCE_FOLDER}/t2tv2.greylist.k29 ${ECTMP_FOLDER}

echo $(date) - Generate genome kmer database..
export ECTMP_FOLDER=${ECRESOURCE_FOLDER}/t2tv2.tmp/
[ ! -d "${ECTMP_FOLDER}" ] && mkdir -p "${ECTMP_FOLDER}"
./pipeline/kmc_3.2.1.sif kmc -k29 -m12 -t12 -ci1 -cs1000000 \
-fm ${ECGENOME_FOLDER}/t2tv2.fasta \
${ECRESOURCE_FOLDER}/t2tv2.k29 ${ECTMP_FOLDER}
echo $(date) - Generate genome minimap2 index..
./pipeline/minimap2_2.24.sif minimap2 -ax map-ont -t 16 \
-d ${ECGENOME_FOLDER}/t2tv2.mmi \
${ECGENOME_FOLDER}/t2tv2.fasta

echo $(date) - Generate MT kmer database..
export ECTMP_FOLDER=${ECRESOURCE_FOLDER}/MT.tmp/
[ ! -d "${ECTMP_FOLDER}" ] && mkdir -p "${ECTMP_FOLDER}"
./pipeline/kmc_3.2.1.sif kmc -k29 -m12 -t12 -ci1 -cs1000000 \
-fm ${ECGENOME_FOLDER}/MT.fasta \
${ECRESOURCE_FOLDER}/MT.k29 ${ECTMP_FOLDER}
echo $(date) - Generate MT minimap2 index..
./pipeline/minimap2_2.24.sif minimap2 -ax map-ont -t 16 \
-d ${ECGENOME_FOLDER}/MT.mmi \
${ECGENOME_FOLDER}/MT.fasta

echo $(date) - End of kmers databases setup.
