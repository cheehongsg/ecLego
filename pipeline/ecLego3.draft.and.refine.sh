#!/bin/bash
set -eu -o pipefail

##
## ecLego3 v3.2404.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##
## ecLego3 workflow : draft, refine and assess the CGs
## This is the top-level driver script
##

export ECLEGO_PIPELINE_PATH=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
export ECLEGO_ROOT=$(dirname "${ECLEGO_PIPELINE_PATH}")
export ECLEGO_RESOURCE_PATH=${ECLEGO_PIPELINE_PATH}/resources

export ECLEGO_RESOURCE_GENOME=${ECLEGO_RESOURCE_PATH}/t2tv2.k29
export ECLEGO_RESOURCE_MT=${ECLEGO_RESOURCE_PATH}/MT.k29

if command -v apptainer &> /dev/null; then
    export APPTAINER_BINDPATH=${ECLEGO_ROOT}
else
    export SINGULARITY_BINDPATH=${ECLEGO_ROOT}
fi

export FASTQ=${1}
export REFGENOME=${2}
export UDFcov2CN=${3:-""}
export UDFSEEDEXT=${4:-""}
export UDFMINOVERLAP=${5:-""}

export SEEDEXT=5000
export MINOVERLAP=$(expr ${SEEDEXT} \* 2)

export KMCTHREADS=20
export KMCMAXMEMORYGB=16
export MINKMERCOUNT=2
export MAXKMERCOUNT=100000
export MINIMAP2THREADS=16
export SAMTOOLSTHREADS=16
export FLYETHREADS=16
export FLYEMINOVERLAP=${MINOVERLAP}

export OUTPUTDIR="${FASTQ}-ecLegoV3"
[ ! -d "${OUTPUTDIR}" ] && mkdir -p "${OUTPUTDIR}"
export PACKAGEDIR="${FASTQ}-results"
[ ! -d "${PACKAGEDIR}" ] && mkdir -p "${PACKAGEDIR}"
export KMCTMPDIR="${OUTPUTDIR}/kmc_tmp"
[ ! -d "${KMCTMPDIR}" ] && mkdir -p "${KMCTMPDIR}"
export OUTPREFIX="${OUTPUTDIR}/$(basename ${FASTQ})"

export RESULTPREFIX=$(basename ${FASTQ})
export RESULTPREFIX=$(echo $RESULTPREFIX | sed 's?.fq.gz$\|.fastq.gz$\|.fq$\|.fastq$??')


# set up the command
export PIGZCMD="singularity exec ${ECLEGO_PIPELINE_PATH}/pigz_2.3.4.sif pigz"
# export PIGZCMD="pigz"
export KMCCMD="${ECLEGO_PIPELINE_PATH}/kmc_3.2.1.sif kmc"
export KMCTOOLSCMD="${ECLEGO_PIPELINE_PATH}/kmc_3.2.1.sif kmc_tools"
export MINIMAP2CMD="${ECLEGO_PIPELINE_PATH}/minimap2_2.24.sif minimap2"
export SAMTOOLSCMD="${ECLEGO_PIPELINE_PATH}/samtools_1.19.sif samtools"
export FLYECMD="${ECLEGO_PIPELINE_PATH}/flye_2.9.1.sif flye"

export ECLGKMCTOOLSCMD="${ECLEGO_PIPELINE_PATH}/kmc_suite_3.2.1.mod.sif kmc_tools"
export THRESHOLDSCMD="perl ${ECLEGO_PIPELINE_PATH}/getThreshold.pl"
export ECLEGO3CMD="perl ${ECLEGO_PIPELINE_PATH}/ecLego3.pl"

export GENOMECOVERAGEBEDCMD="${ECLEGO_PIPELINE_PATH}/bedtools_2.31.0.sif genomeCoverageBed"
export SORTBEDCMD="${ECLEGO_PIPELINE_PATH}/bedtools_2.31.0.sif sortBed"
export BG2BWCMD="${ECLEGO_PIPELINE_PATH}/ucsc-bedgraphtobigwig_445.sif bedGraphToBigWig"


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
export MINAMPLICATIONCN=${cov5CN}
echo $(date) Thresholding COMPLETED

#
if [[ -n "${UDFcov2CN}" ]]; then
    cov2cn_value=$(echo "${UDFcov2CN}" | bc)
    if [[ $cov2cn_value -gt 0 ]]; then
        echo "# Adjusting Cov2CN ${cov2CN} to user-defined value ${cov2cn_value}"
        cov2CN="${cov2cn_value}"
        covcn_value=$(expr ${cov2cn_value} / 2)
        echo "# Adjusting CovCN ${covCN} to user-defined value ${covcn_value}"
        covCN="${covcn_value}"
        echo "# Not adjusting the remaining coverage variables."
    else
        echo "#####"
        echo "# WARNING: Ignoring unusable coverage threshold provided ${UDFcov2CN}"
        echo "#####"
    fi
fi

if [[ -n "${UDFSEEDEXT}" ]]; then
    seedext_value=$(echo "${UDFSEEDEXT}" | bc)
    if [[ $seedext_value -gt 0 ]]; then
        echo "# Adjusting seed extension ${SEEDEXT} to user-defined value ${seedext_value}"
        SEEDEXT="${seedext_value}"
    else
        echo "#####"
        echo "# WARNING: Ignoring unusable seed extension provided ${UDFSEEDEXT}"
        echo "#####"
    fi
fi

if [[ -n "${UDFMINOVERLAP}" ]]; then
    minoverlap_value=$(echo "${UDFMINOVERLAP}" | bc)
    if [[ $minoverlap_value -gt 0 ]]; then
        MINOVERLAP=${minoverlap_value}
        echo "# Adjusting minimum overlap ${FLYEMINOVERLAP} accordingly to ${MINOVERLAP}"
        FLYEMINOVERLAP=${MINOVERLAP}
    else
        echo "#####"
        echo "# WARNING: Ignoring unusable minimum overlap provided ${UDFMINOVERLAP}"
        echo "#####"
    fi
fi

# report final values
echo "# Effective values ---"
echo "# cov2CN=${cov2CN}"
echo "# covCN=${covCN}"
echo "# SEEDEXT=${SEEDEXT}"
echo "# MINOVERLAP=${MINOVERLAP}"
echo "# FLYEMINOVERLAP=${FLYEMINOVERLAP}"
echo "# ---"

##### TODO: package the custom kmc_tools
echo $(date) Amplified reads selection STARTED
time (${ECLGKMCTOOLSCMD} -v -t12 ampreads \
${OUTPREFIX} -cn${covCN} \
nokmcdb \
${FASTQ} \
-fq -ci2 ${OUTPREFIX}.readamp.cn_${covCN}.fastq \
1>${OUTPREFIX}.readamp.cn_${covCN}.fastq.stdout \
2>${OUTPREFIX}.readamp.cn_${covCN}.fastq.stderr)
echo $(date) Amplified reads selection ENDED
# 

#
wc -l ${OUTPREFIX}.readamp.cn_${covCN}.fastq
#


# NOTE: while waiting for the long-process of 1st filter on readamp.cn_${covCN}.fastq
#       we proceed with downstream testing!!!
export KMCMAXMEMORYGB=16
export KMCTHREADS=12
export MAXKMERCOUNT=100000
# generate the read-depth-count (RDC) of the candidate set
echo $(date) Amplified reads kmerdb generation STARTED
time (/net/nwgc/vol1/techdev/GBMTW/KMC-3.2.1/bin/kmc \
-k29 -m${KMCMAXMEMORYGB} -t${KMCTHREADS} -cs${MAXKMERCOUNT} \
-fq ${OUTPREFIX}.readamp.cn_${covCN}.fastq \
${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc \
${KMCTMPDIR} 2>&1 | tee ${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.log)
echo $(date) Amplified reads kmerdb generation ENDED
#

echo $(date) Baseline amplified reads exclusion STARTED
time (${ECLGKMCTOOLSCMD} -v -t12 ampreads \
${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc -cn${covCN} \
${ECLEGO_RESOURCE_GENOME} \
${OUTPREFIX}.readamp.cn_${covCN}.fastq \
-fq -ci3 ${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.t2tv2.cn3.fastq \
1>${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.t2tv2.cn3.fastq.stdout \
2>${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.t2tv2.cn3.fastq.stderr)
echo $(date) Baseline amplified reads exclusion ENDED
#


###
# iteration #1
###

# this will over-filtered, but most of the noise are gone!
# let's now recover those that have been over-filtered
# i.e. we pick up reads that overlap the seeds kmcdb
#
# build the kmcdb for quick picking
#
echo $(date) Seeds kmerdb generation STARTED
time (/net/nwgc/vol1/techdev/GBMTW/KMC-3.2.1/bin/kmc \
-k29 -m${KMCMAXMEMORYGB} -t${KMCTHREADS} -cs${MAXKMERCOUNT} \
-fq ${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.t2tv2.cn3.fastq \
${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.t2tv2.cn3.fastq.rdc \
${KMCTMPDIR} 2>&1 | tee ${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.t2tv2.cn3.fastq.rdc.log)
echo $(date) Seeds kmerdb generation ENDED



echo $(date) Seeds extension STARTED
time (${ECLGKMCTOOLSCMD} -t12 filter \
${OUTPREFIX}.readamp.cn_${covCN}.fastq.rdc.t2tv2.cn3.fastq.rdc \
-ci1 ${OUTPREFIX}.readamp.cn_${covCN}.fastq \
-fq -ci${SEEDEXT} ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.fastq)
echo $(date) Seeds extension ENDED

#
wc -l ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.fastq
# 



echo $(date) Greylist purging STARTED
export GENOMEDIR=$(dirname ${REFGENOME})
time (${ECLGKMCTOOLSCMD} -t12 filter \
${GENOMEDIR}/../resources/t2tv2.greylist.k29 \
-ci2 ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.fastq \
-fq -cx0.409 ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq)
#
echo $(date) Greylist purging ENDED
wc -l ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq




#####
# finally generate the MT and noMT

# MT only
echo $(date) MT candidate fastq generation STARTED
time (${ECLGKMCTOOLSCMD} -t12 filter \
${ECLEGO_RESOURCE_MT} -ci1 \
${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq \
-fq -ci${SEEDEXT} ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq.MT.fastq)
#
#
wc -l ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq.MT.fastq
#
ln ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq.MT.fastq \
${OUTPREFIX}.MT.ci${SEEDEXT}.fastq
time (${PIGZCMD} -f -p 4 ${OUTPREFIX}.MT.ci${SEEDEXT}.fastq)
echo $(date) MT candidate fastq generation ENDED
#

# nonMT only
echo $(date) Amplified candidate fastq generation STARTED
time (${ECLGKMCTOOLSCMD} -t12 filter \
${ECLEGO_RESOURCE_MT} -ci1 \
${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq \
-fq -ci0 -cx4999 ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq.nonMT.fastq)
#
wc -l ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq.nonMT.fastq
#
ln ${OUTPREFIX}.readamp.cn_${covCN}.fastq.iter2.ci${SEEDEXT}.nogreylist.fastq.nonMT.fastq \
${OUTPREFIX}.candidates.ci${SEEDEXT}.fastq
time (${PIGZCMD} -f -p 4 ${OUTPREFIX}.candidates.ci${SEEDEXT}.fastq)
echo $(date) Amplified candidate fastq generation ENDED
#



# perform de novo on MT
# --extra-params max_coverage_drop_rate=2,max_extensions_drop_rate=2 \
echo $(date) de Novo MT assembly STARTED
${FLYECMD} --threads ${FLYETHREADS} --meta \
--nano-hq ${OUTPREFIX}.MT.ci${SEEDEXT}.fastq.gz \
--out-dir ${OUTPREFIX}_MT_assembly.ci${SEEDEXT} \
--extra-params max_coverage_drop_rate=2,max_extensions_drop_rate=2 \
-m ${FLYEMINOVERLAP} &
wait
echo $(date) de Novo MT assembly COMPLETED
#




# perform de novo on candidates
# --extra-params max_coverage_drop_rate=2,max_extensions_drop_rate=2 \
echo $(date) de Novo candidates assembly STARTED
${FLYECMD} --threads ${FLYETHREADS} --meta \
--nano-hq ${OUTPREFIX}.candidates.ci${SEEDEXT}.fastq.gz \
--out-dir ${OUTPREFIX}_candidate_assembly.ci${SEEDEXT} \
--extra-params max_coverage_drop_rate=2,max_extensions_drop_rate=2 \
-m ${FLYEMINOVERLAP} &
wait
echo $(date) de Novo candidates assembly COMPLETED
#




# optional: perform mapping as checking {graph_before_rr, graph_after_rr, assembly}
export GRAPHFASTA=${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/20-repeat/graph_before_rr.fasta ; \
[ ! -z "${REFGENOME}" ] && [ -f "${REFGENOME}" ] && [ -f "${GRAPHFASTA}" ] && { 
# graph_before_rr
echo $(date) Assembly graph_before_rr reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${REFGENOME}" \
${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/20-repeat/graph_before_rr.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${OUTPREFIX}.candidates.ci${SEEDEXT}.graph_before_rr.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${OUTPREFIX}.candidates.ci${SEEDEXT}.graph_before_rr.bam ; \
export DSNAME=${OUTPREFIX}.candidates.ci${SEEDEXT}.graph_before_rr.bam;
export CNSCALE=$(awk -v DICNSCALE=${cov2CN} 'BEGIN{print 2/DICNSCALE}') ; \
echo $(date) Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} \
| ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo $(date) Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg ${REFGENOME}.fai ${DSNAME}.2cn_${cov2CN}.bw && rm ${DSNAME}.bg ; \
echo $(date) Done. ) &
wait
echo $(date) Assembly graph_before_rr reference mapping COMPLETED

# filtered+raw summary
echo $(date) ecDNAs summarizing STARTED
${ECLEGO3CMD} \
sievegraph \
--gv ${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/20-repeat/graph_before_rr.gv \
--diploidcov ${cov2CN} \
--minlenkb 25 \
--oprefix ${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/${RESULTPREFIX} \
--bam ${OUTPREFIX}.candidates.ci${SEEDEXT}.graph_before_rr.bam
echo $(date) ecDNAs summarizing ENDED

}
#






# optional: perform mapping as checking {graph_before_rr, graph_after_rr, assembly}
export GRAPHFASTA=${OUTPREFIX}_MT_assembly.ci${SEEDEXT}/20-repeat/graph_before_rr.fasta ; \
[ ! -z "${REFGENOME}" ] && [ -f "${REFGENOME}" ] && [ -f "${GRAPHFASTA}" ] && { 
# graph_before_rr
echo $(date) Assembly graph_before_rr reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${REFGENOME}" \
${OUTPREFIX}_MT_assembly.ci${SEEDEXT}/20-repeat/graph_before_rr.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@4 -o ${OUTPREFIX}.MT.ci${SEEDEXT}.graph_before_rr.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${OUTPREFIX}.MT.ci${SEEDEXT}.graph_before_rr.bam ; \
export DSNAME=${OUTPREFIX}.MT.ci${SEEDEXT}.graph_before_rr.bam;
export CNSCALE=$(awk -v DICNSCALE=${cov2CN} 'BEGIN{print 2/DICNSCALE}') ; \
echo $(date) Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} \
| ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo $(date) Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg ${REFGENOME}.fai ${DSNAME}.2cn_${cov2CN}.bw && rm ${DSNAME}.bg ; \
echo $(date) Done. ) &
wait
echo $(date) Assembly graph_before_rr reference mapping COMPLETED

# summary
echo $(date) MT summarizing STARTED
${ECLEGO3CMD} \
sievegraph \
--gv ${OUTPREFIX}_MT_assembly.ci${SEEDEXT}/20-repeat/graph_before_rr.gv \
--diploidcov ${cov2CN} \
--oprefix ${OUTPREFIX}_MT_assembly.ci${SEEDEXT}/${RESULTPREFIX} \
--bam ${OUTPREFIX}.MT.ci${SEEDEXT}.graph_before_rr.bam
echo $(date) MT summarizing ENDED

}
#

echo $(date) Generate refinement script STARTED
${ECLEGO3CMD} \
initsubgraph \
--summary ${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/${RESULTPREFIX}.cn${cov2CN}.disjointcyclic.gv.overview.xls \
--genome ${REFGENOME} \
--scriptPath ${ECLEGO_PIPELINE_PATH} \
--diploidcov ${cov2CN} \
--seedext ${SEEDEXT} \
--minoverlap ${FLYEMINOVERLAP} \
--step refine \
> ecLego3.step.refine.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
chmod a+x ecLego3.step.refine.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo $(date) Generate refinement script ENDED


echo $(date) Generate refinement assessment script STARTED
${ECLEGO3CMD} \
initsubgraph \
--summary ${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/${RESULTPREFIX}.cn${cov2CN}.disjointcyclic.gv.overview.xls \
--genome ${REFGENOME} \
--scriptPath ${ECLEGO_PIPELINE_PATH} \
--diploidcov ${cov2CN} \
--seedext ${SEEDEXT} \
--minoverlap ${FLYEMINOVERLAP} \
--step assess \
> ecLego3.step.assess.refined.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
chmod a+x ecLego3.step.assess.refined.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo $(date) Generate refinement assessment script ENDED


#####
# Draft assembly refinement handling
#####


echo 
echo '#####'
echo '#' $(date) Refinement of draft assembly STARTED
echo '#####'
./ecLego3.step.refine.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo '#####'
echo '#' $(date) Refinement of draft assembly ENDED
echo '#####'


#####
# SVs handling
#####


echo 
echo '#####'
echo '#' $(date) Assessment of refined assembly STARTED
echo '#####'
./ecLego3.step.assess.refined.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo '#####'
echo '#' $(date) Assessment of refined assembly ENDED
echo '#####'


#####
# subpopulation handling
#####


echo $(date) Subpopulations script STARTED
${ECLEGO3CMD} \
initsubpop \
--summary ${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/${RESULTPREFIX}.cn${cov2CN}.disjointcyclic.gv.overview.xls \
--genome ${REFGENOME} \
--scriptPath ${ECLEGO_PIPELINE_PATH} \
--diploidcov ${cov2CN} \
--seedext ${SEEDEXT} \
--minoverlap ${FLYEMINOVERLAP} \
> ecLego3.step.generate.subpopulations.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
chmod a+x ecLego3.step.generate.subpopulations.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo $(date) Subpopulations script ENDED


echo 
echo '#####'
echo '#' $(date) Subpopulations data files generation STARTED
echo '#####'
./ecLego3.step.generate.subpopulations.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo '#####'
echo '#' $(date) Subpopulations data files generation ENDED
echo '#####'


#####
# harvest the results for downstream analysis and manual inspection
#####


echo $(date) Harvest results script STARTED
${ECLEGO3CMD} \
harvestresult \
--summary ${OUTPREFIX}_candidate_assembly.ci${SEEDEXT}/${RESULTPREFIX}.cn${cov2CN}.disjointcyclic.gv.overview.xls \
--genome ${REFGENOME} \
--scriptPath ${ECLEGO_PIPELINE_PATH} \
--diploidcov ${cov2CN} \
--seedext ${SEEDEXT} \
--minoverlap ${FLYEMINOVERLAP} \
--resultPath ${PACKAGEDIR} \
> ecLego3.step.harvest.results.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
chmod a+x ecLego3.step.harvest.results.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo $(date) Harvest results scriptENDED


echo 
echo '#####'
echo '#' $(date) Harvest results STARTED
echo '#####'
./ecLego3.step.harvest.results.se${SEEDEXT}.dicn${cov2CN}.ovlp${FLYEMINOVERLAP}.sh
echo '#####'
echo '#' $(date) Harvest results ENDED
echo '#####'

