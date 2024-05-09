#!/bin/bash
set -eu -o pipefail

##
## ecLego2 v2.04.01
## Copyright (C) 2019-2023 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

echo $(date) - Downloading sif files to pipeline subfolder..
curl -o pipeline/samtools_v1.15.1.sif https://depot.galaxyproject.org/singularity/samtools%3A1.15.1--h1170115_0 -C -
curl -o pipeline/seqtk_1.3.sif https://depot.galaxyproject.org/singularity/seqtk%3A1.3--h7132678_4 -C -
curl -o pipeline/ucsc-bigbedtobed_377.sif https://depot.galaxyproject.org/singularity/ucsc-bigbedtobed%3A377--ha8a8165_3 -C -
curl -o pipeline/bedtools_2.31.0.sif https://depot.galaxyproject.org/singularity/bedtools%3A2.31.0--hf5e1c6e_2 -C -
curl -o pipeline/kmc_3.2.1.sif https://depot.galaxyproject.org/singularity/kmc%3A3.2.1--hf1761c0_2 -C -
curl -o pipeline/flye_2.9.1.sif https://depot.galaxyproject.org/singularity/flye%3A2.9.1--py39h6935b12_0 -C -
curl -o pipeline/pigz_2.3.4.sif https://depot.galaxyproject.org/singularity/pigz%3A2.3.4 -C -
curl -o pipeline/minimap2_2.24.sif https://depot.galaxyproject.org/singularity/minimap2%3A2.24--h7132678_1 -C -
curl -o pipeline/igvtools_2.14.1.sif https://depot.galaxyproject.org/singularity/igvtools%3A2.14.1--hdfd78af_0 -C -
echo $(date) - End of download.

echo $(date) - Set up sif execution permission..
chmod u+x pipeline/*.sif
echo $(date) - End of permission set up.
