#!/bin/bash
set -eu -o pipefail

##
## ecLego2 v2.04.01
## Copyright (C) 2019-2023 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##
## example:
##   assume the folder structure as such:
##   <path_to_ecLego2_root_folder>/
##       pipeline/
##           genomes/
##               t2tv2.fasta
##               t2tv2.fasta.fai
##           resources/
##               t2tv2.greylist.k29.kmc_pre
##               t2tv2.greylist.k29.kmc_suf
##       data/
##           A01.sup/
##               A01.sup.fastq.gz
##
##   export SINGULARITY_BINDPATH=<path_to_ecLego2_root_folder>
##   cd <path_to_ecLego2_root_folder>
##   pipeline/ecLegov2_s01.sh data/A01.sup/A01.sup.fastq.gz pipeline/genomes/t2tv2.fasta 2>&1 | tee data/A01.sup/A01.sup.fastq.gz.ecLegov2.log
##


export FASTQ=${1}
export REFGENOME=${2}
export KMCTHREADS=${SLURM_CPUS_PER_TASK}
export KMCMAXMEMORYGB=16
export MINKMERCOUNT=2
export MAXKMERCOUNT=100000
export MINIMAP2THREADS=${SLURM_CPUS_PER_TASK}
export SAMTOOLSTHREADS=${SLURM_CPUS_PER_TASK}
export FLYETHREADS=${SLURM_CPUS_PER_TASK}
export FLYEMINOVERLAP=10000

# check for empty parameters!
[ -z "${FASTQ}" ] && { echo "FASTQ is empty" ; exit 1; }
[ ! -z "${REFGENOME}" ] && [ ! -f "${REFGENOME}" ] && { echo "Reference genome is empty" ; exit 1; }

# owner and group access; block all others
umask 0007

echo $(date) ecLegov2 STARTED
# report parameters
echo FASTQ=${FASTQ}
echo KMCTHREADS=${KMCTHREADS}
echo MINKMERCOUNT=${MINKMERCOUNT}
echo MAXKMERCOUNT=${MAXKMERCOUNT}
echo SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}

# check parameters
[ ! -f "${FASTQ}" ] && { echo "${FASTQ} does not exist" ; exit 2 ; }

# create output directory
export OUTPUTDIR="${FASTQ}-output"
[ ! -d "${OUTPUTDIR}" ] && mkdir -p "${OUTPUTDIR}"
export OUTPREFIX="${OUTPUTDIR}/$(basename ${FASTQ})"
export RESULTPREFIX=$(basename ${FASTQ})
export RESULTPREFIX=${RESULTPREFIX/%.fastq.gz/}

# set up the command
export KMCCMD="./pipeline/kmc_3.2.1.sif kmc"
export KMCTOOLSCMD="./pipeline/kmc_3.2.1.sif kmc_tools"
export THRESHOLDSCMD="perl ./pipeline/getThreshold.pl"
export MINIMAP2CMD="./pipeline/minimap2_2.24.sif minimap2"
export SAMTOOLSCMD="./pipeline/samtools_v1.15.1.sif samtools"
export IGVTOOLSCMD="./pipeline/igvtools_2.14.1.sif igvtools"
export FLYECMD="./pipeline/flye_2.9.1.sif flye"

[ ! -f "$(echo $KMCCMD | cut -f 1 -d\  )" ] && { echo "$(echo $KMCCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $KMCTOOLSCMD | cut -f 1 -d\  )" ] && { echo "$(echo $KMCTOOLSCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $MINIMAP2CMD | cut -f 1 -d\  )" ] && { echo "$(echo $MINIMAP2CMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $SAMTOOLSCMD | cut -f 1 -d\  )" ] && { echo "$(echo $SAMTOOLSCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $IGVTOOLSCMD | cut -f 1 -d\  )" ] && { echo "$(echo $IGVTOOLSCMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }
[ ! -f "$(echo $FLYECMD | cut -f 1 -d\  )" ] && { echo "$(echo $FLYECMD | cut -f 1 -d\  ) does not exist, SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH}" ; exit 2 ; }

# construct k-mer spectrum
echo $(date) K-mer construction STARTED
${KMCCMD} -k29 -m${KMCMAXMEMORYGB} -t${KMCTHREADS} -cs${MAXKMERCOUNT} \
-fq ${FASTQ} \
${OUTPREFIX} \
${OUTPUTDIR}/ &
wait
echo $(date) K-mer construction COMPLETED

# report k-mer spectrum
echo $(date) K-mer report STARTED
${KMCTOOLSCMD} transform \
${OUTPREFIX} \
histogram ${OUTPREFIX}.hist -ci${MINKMERCOUNT} -cx${MAXKMERCOUNT} &
wait
echo $(date) K-mer report COMPLETED

# set up the thresholds
echo $(date) Thresholding STARTED
${THRESHOLDSCMD} ${OUTPREFIX}.hist
export $(${THRESHOLDSCMD} ${OUTPREFIX}.hist | grep -v '^#' | xargs)
echo $(date) Thresholding COMPLETED

# filter1: CNA reads
echo $(date) CNA filtering STARTED
${KMCTOOLSCMD} filter \
${OUTPREFIX} -ci${cov5CN} ${FASTQ} \
-fq -ci0.41 ${OUTPREFIX}.filter1.fastq &
wait
echo $(date) CNA filtering COMPLETED

# filter2: chromosomal reads
echo $(date) Chromosomal filtering STARTED
${KMCTOOLSCMD} filter \
${OUTPREFIX} -ci${covCN} -cx$(expr ${cov5CN} - 1) ${OUTPREFIX}.filter1.fastq \
-fq -cx${errorRate} ${OUTPREFIX}.filter2.fastq &
wait
echo $(date) Chromosomal filtering COMPLETED

# compress filter1
echo $(date) Compressing filter1 fastq STARTED
gzip ${OUTPREFIX}.filter1.fastq &
wait
echo $(date) Compressing filter1 fastq COMPLETED

# compress filter2
echo $(date) Compressing filter2 fastq STARTED
gzip ${OUTPREFIX}.filter2.fastq &
wait
echo $(date) Compressing filter2 fastq COMPLETED

# filter3: greylist
echo $(date) Final filtering STARTED
${KMCTOOLSCMD} filter \
pipeline/resources/t2tv2.greylist.k29 -ci2 ${OUTPREFIX}.filter2.fastq.gz \
-fq -cx0.409 ${OUTPREFIX}.filter3.fastq &
wait
echo $(date) Final filtering COMPLETED

# compress filter3
echo $(date) Compressing filter3 fastq STARTED
gzip ${OUTPREFIX}.filter3.fastq &
wait
echo $(date) Compressing filter3 fastq COMPLETED

# second distillation round
echo $(date) Distillation k-mer construction STARTED
${KMCCMD} -k29 -m${KMCMAXMEMORYGB} -t${KMCTHREADS} -ci1 -cs${MAXKMERCOUNT} \
-fq ${OUTPREFIX}.filter3.fastq.gz \
${OUTPREFIX}.2ndpass \
${OUTPUTDIR}/ &
wait
echo $(date) Distillation k-mer construction COMPLETED

# expected
echo $(date) Distillation baselining STARTED
${KMCTOOLSCMD} simple \
${OUTPREFIX}.2ndpass -ci${cov5CN} \
pipeline/resources/t2tv2.k29 \
intersect \
${OUTPREFIX}.2ndpass.baseline -ocright &
wait
echo $(date) Distillation baselining COMPLETED
#
echo $(date) Writing distillation baseline STARTED
${KMCTOOLSCMD} transform \
${OUTPREFIX}.2ndpass.baseline \
dump ${OUTPREFIX}.2ndpass.baseline.txt &
wait
echo $(date) Writing distillation baseline COMPLETED

# empirical
echo $(date) Distillation sampling STARTED
${KMCTOOLSCMD} simple \
${OUTPREFIX}.2ndpass -ci${cov5CN} \
pipeline/resources/t2tv2.k29 \
intersect \
${OUTPREFIX}.2ndpass.sample -ocleft &
wait
echo $(date) Distillation sampling COMPLETED
#
echo $(date) Writing distillation sample STARTED
${KMCTOOLSCMD} transform \
${OUTPREFIX}.2ndpass.sample \
dump ${OUTPREFIX}.2ndpass.sample.txt &
wait
echo $(date) Writing distillation sample COMPLETED

# empirical
echo $(date) Distillation unique STARTED
${KMCTOOLSCMD} simple \
${OUTPREFIX}.2ndpass -ci${cov5CN} \
pipeline/resources/t2tv2.k29 \
kmers_subtract \
${OUTPREFIX}.2ndpass.unique &
wait
echo $(date) Distillation unique COMPLETED
#
echo $(date) Writing distillation unique STARTED
${KMCTOOLSCMD} transform \
${OUTPREFIX}.2ndpass.unique \
dump ${OUTPREFIX}.2ndpass.unique.txt &
wait
echo $(date) Writing distillation unique COMPLETED

# collate
(paste ${OUTPREFIX}.2ndpass.sample.txt ${OUTPREFIX}.2ndpass.baseline.txt \
| awk -v CI=${cov2CN} -v XC=${cov6CN} -v FS="\t" -v OFS="\t" '{if (($2-($4*CI))>=XC) { print ">"$1"\n"$1;}}'; \
awk -v FS="\t" -v OFS="\t" '{ print ">"$1"\n"$1;}' ${OUTPREFIX}.2ndpass.unique.txt ) \
> ${OUTPREFIX}.2ndpass.fasta

# compress files
echo $(date) Compressing distillation file\# 1 STARTED
gzip ${OUTPREFIX}.2ndpass.sample.txt &
wait
echo $(date) Compressing distillation file\# 1 COMPLETED
echo $(date) Compressing distillation file\# 2 STARTED
gzip ${OUTPREFIX}.2ndpass.baseline.txt &
wait
echo $(date) Compressing distillation file\# 2 COMPLETED
echo $(date) Compressing distillation file\# 3 STARTED
gzip ${OUTPREFIX}.2ndpass.unique.txt &
wait
echo $(date) Compressing distillation file\# 3 COMPLETED

# second distillation round
echo $(date) Distilled k-mer construction STARTED
${KMCCMD} -k29 -m${KMCMAXMEMORYGB} -t${KMCTHREADS} -ci1 -cs${MAXKMERCOUNT} \
-fa ${OUTPREFIX}.2ndpass.fasta \
${OUTPREFIX}.distilled \
${OUTPUTDIR}/ &
wait
echo $(date) Distilled k-mer construction COMPLETED

echo $(date) Compressing distillation file\# 4 STARTED
gzip ${OUTPREFIX}.2ndpass.fasta &
wait
echo $(date) Compressing distillation file\# 4 COMPLETED

# write final set of reads
echo $(date) Finalizing distilled reads STARTED
${KMCTOOLSCMD} filter \
${OUTPREFIX}.distilled -ci1 \
${OUTPREFIX}.filter3.fastq.gz \
-fq -ci0.41 ${OUTPREFIX}.candidates.fastq &
wait
echo $(date) Finalizing distilled reads COMPLETED

# compress distilled reads
echo $(date) Compressing distilled reads fastq STARTED
gzip ${OUTPREFIX}.candidates.fastq &
wait
echo $(date) Compressing distilled reads fastq COMPLETED

# optional: perform mapping as checking
[ ! -z "${REFGENOME}" ] && [ -f "${REFGENOME}" ] && { 
echo $(date) Reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} --secondary=no \
"${REFGENOME}" \
${OUTPREFIX}.candidates.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${OUTPREFIX}.candidates.bam
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${OUTPREFIX}.candidates.bam ; \
${IGVTOOLSCMD} count -w 1 -z 7 \
${OUTPREFIX}.candidates.bam \
${OUTPREFIX}.candidates.tdf \
"${REFGENOME}" ) &
wait
echo $(date) Reference mapping COMPLETED
}

# perform de novo
echo $(date) de Novo assembly STARTED
${FLYECMD} --threads ${FLYETHREADS} --meta \
--nano-hq ${OUTPREFIX}.candidates.fastq.gz \
--out-dir ${OUTPREFIX}_candidate_assembly \
-m ${FLYEMINOVERLAP} &
wait
echo $(date) de Novo assembly COMPLETED

# optional: perform mapping as checking {graph_before_rr, graph_after_rr, assembly}
[ ! -z "${REFGENOME}" ] && [ -f "${REFGENOME}" ] && { 
# graph_before_rr
echo $(date) Assembly graph_before_rr reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${REFGENOME}" \
${OUTPREFIX}_candidate_assembly/20-repeat/graph_before_rr.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${OUTPREFIX}.candidates.graph_before_rr.bam
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${OUTPREFIX}.candidates.graph_before_rr.bam ; \
${IGVTOOLSCMD} count -w 1 -z 7 \
${OUTPREFIX}.candidates.graph_before_rr.bam \
${OUTPREFIX}.candidates.graph_before_rr.tdf \
"${REFGENOME}" ) &
wait
echo $(date) Assembly graph_before_rr reference mapping COMPLETED

# summary
echo $(date) ecDNAs summarizing STARTED
perl ./pipeline/ecLegov2.pl \
sievegraph \
--gv ${OUTPREFIX}_candidate_assembly/20-repeat/graph_before_rr.gv \
--diploidcov ${cov2CN} \
--oprefix ${OUTPREFIX}_candidate_assembly/${RESULTPREFIX} \
--bam ${OUTPREFIX}.candidates.graph_before_rr.bam
echo $(date) ecDNAs summarizing ENDED

[ -f "${OUTPREFIX}_candidate_assembly/20-repeat/graph_before_rr.fasta" ] && { 
# graph_before_rr as reference
echo $(date) Reads to graph_before_rr mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
${OUTPREFIX}_candidate_assembly/20-repeat/graph_before_rr.fasta \
${OUTPREFIX}.candidates.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${OUTPREFIX}.candidates.graph_before_rr.asref.bam
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${OUTPREFIX}.candidates.graph_before_rr.asref.bam ; \
${IGVTOOLSCMD} count -w 1 -z 7 \
${OUTPREFIX}.candidates.graph_before_rr.asref.bam \
${OUTPREFIX}.candidates.graph_before_rr.asref.tdf \
"${OUTPREFIX}_candidate_assembly/20-repeat/graph_before_rr.fasta" ) &
wait
echo $(date) Reads to graph_before_rr COMPLETED
}

# graph_after_rr
echo $(date) Assembly graph_after_rr reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${REFGENOME}" \
${OUTPREFIX}_candidate_assembly/20-repeat/repeat_graph_edges.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${OUTPREFIX}.candidates.graph_after_rr.bam
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${OUTPREFIX}.candidates.graph_after_rr.bam ; \
${IGVTOOLSCMD} count -w 1 -z 7 \
${OUTPREFIX}.candidates.graph_after_rr.bam \
${OUTPREFIX}.candidates.graph_after_rr.tdf \
"${REFGENOME}" ) &
wait
echo $(date) Assembly graph_after_rr reference mapping COMPLETED
# graph_after_rr
echo $(date) Assembly reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${REFGENOME}" \
${OUTPREFIX}_candidate_assembly/assembly.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${OUTPREFIX}.candidates.assembly.bam
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${OUTPREFIX}.candidates.assembly.bam ; \
${IGVTOOLSCMD} count -w 1 -z 7 \
${OUTPREFIX}.candidates.assembly.bam \
${OUTPREFIX}.candidates.assembly.tdf \
"${REFGENOME}" ) &
wait
echo $(date) Assembly reference mapping COMPLETED
}

#
echo $(date) ecLegov2 ENDED
