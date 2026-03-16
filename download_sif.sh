#!/bin/bash
set -eu -o pipefail

##
## ecLego3 v3.2404.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

echo $(date) - Downloading sif files to pipeline subfolder..
# for database setup
curl -# -L -o pipeline/samtools_1.19.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/samtools%3A1.19--h50ea8bc_0 -C -
curl -# -L -o pipeline/ucsc-bigbedtobed_377.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/ucsc-bigbedtobed%3A377--ha8a8165_3 -C -
curl -# -L -o pipeline/bedtools_2.31.0.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/bedtools%3A2.31.0--hf5e1c6e_2 -C -
curl -# -L -o pipeline/seqtk_1.3.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/seqtk%3A1.3--h7132678_4 -C -
curl -# -L -o pipeline/kmc_3.2.1.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/kmc%3A3.2.1--hf1761c0_2 -C -
# for ecLego3
### sif files already downloaded are prefixed with '###'
### kmc_3.2.1.sif
### samtools_1.19.sif
### seqtk_1.3.sif
### bedtools_2.31.0.sif
curl -# -L -o pipeline/pigz_2.3.4.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/pigz%3A2.3.4 -C -
curl -# -L -o pipeline/flye_2.9.1.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/flye%3A2.9.1--py39h6935b12_0 -C -
curl -# -L -o pipeline/minimap2_2.24.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/minimap2%3A2.24--h7132678_1 -C -
curl -# -L -o pipeline/ucsc-bedgraphtobigwig_445.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig%3A445--h954228d_0 -C -
curl -# -L -o pipeline/bcftools_1.19.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_0 -C -
curl -# -L -o pipeline/sniffles_2.3.3.sif -w "%{filename_effective}\n" https://depot.galaxyproject.org/singularity/sniffles%3A2.3.3--pyhdfd78af_0 -C -
curl -# -L -o pipeline/kmc_suite_3.2.1.mod.sif -w "%{filename_effective}\n" 'https://www.dropbox.com/scl/fi/qfb2n1uwkjshtx39mpjx7/kmc_suite_3.2.1.mod.sif?rlkey=jtdwfzmrb6yhvm4h7jq7ua2la&st=ymckleyf&dl=1' -C -
echo $(date) - End of download.

echo $(date) - Set up sif execution permission..
chmod u+x pipeline/*.sif
echo $(date) - End of permission set up.
