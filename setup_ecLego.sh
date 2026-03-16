#!/bin/bash
set -eu -o pipefail

##
## ecLego3 v3.2404.01
## Copyright (C) 2019-2026 Chee-Hong WONG (Dr Chia-Lin Wei Laboratory)
##

echo '#####'
echo '#' $(date) - ecLego system setup - STARTED
echo '#####'

# download singularity images
./download_sif.sh 2>&1 | tee setup.step1.download_sif.run.log

# download t2tv2 genome reference
./download_hsa_t2tv2.sh 2>&1 | tee setup.step2.download_hsa_t2tv2.run.log

# setup kmers databases
./setup_hsa_t2tv2_kmerdbs.sh 2>&1 | tee setup.step3.setup_hsa_t2tv2_kmerdbs.run.log

echo '#####'
echo '#' $(date) - ecLego system setup - ENDED
echo '#####'
