# ecLego3 - circular genomes asssembly from long-read data

To advance ecDNA structural characterization, we devised ecLego3, a customized circular genome assembler workflow, to comprehensively reconstruct ecDNA molecular structures and dissect their diversity from long-read whole genome sequence (lrWGS) data.

## ecLego3 quickstart guide

For detailed documentation and advanced configuration, please visit the [Wiki documentation](../../wiki).


### Installation

Clone the repository into your project directory (e.g., ```/my_project```) and run the setup script:

```bash
# Define your project directory
export PROJECTDIR=/my_project
cd ${PROJECTDIR}

# clone the repository
git clone --branch 'v3.2403.01' --depth 1 https://github.com/cheehongsg/ecLego.git

# Make scripts executable and initialize environment
cd ${PROJECTDIR}/ecLego && chmod u+x ./*.sh
./setup_ecLego.sh
```


### Long-read data setup

By default, ecLego expects long-read FASTQ files to be nested in sample-specific subdirectories under ```/my_project/data```. Following this convention avoids the need for custom container mount configurations.


#### Directory Template:

```
/my_project/data/
└── [Sample_ID]/
    └── [Sample_ID].fastq.gz
```


### Running ecLego3

**Prerequisite:** Ensure Singularity or Apptainer is installed.


#### 1. Generate Read Metrics

First, assess the read length distribution to determine the optimal parameters for your sample.

```bash
export SAMPLEID=my_sample
cd ${PROJECTDIR}/data/${SAMPLEID}

# Calculate FASTQ metrics
${PROJECTDIR}/ecLego/pipeline/getFastqMetrics.pl \
  ${SAMPLEID}.fastq.gz \
  | tee ${SAMPLEID}.fq.gz.metrics
```


#### 2. Execute ecLego3 Pipeline

Apply your derived metrics to initiate the processing script, tailoring parameters as recommended for your sample.

```bash
# Run drafting and refinement
cd ${PROJECTDIR}/data/${SAMPLEID}
${PROJECTDIR}/ecLego/pipeline/ecLego3.draft.and.refine.sh \
  ${SAMPLEID}.fastq.gz \
  ${PROJECTDIR}/ecLego/pipeline/genomes/t2tv2.fasta \
  2>&1 | tee ${SAMPLEID}.ecLego3.run.log
```

**Quick Tip:** Monitoring the ```.log``` file in real-time allows you to catch potential configuration errors early in the drafting phase.


### Parameter Reference

| Position | Parameter | Description | Example Value |
|:---:|---|---|---|
| 1 | Input FASTQ | Path to the long-read sequencing data (compressed). | sample.fastq.gz |
| 2 | Reference FASTA | Path to the reference genome assembly. For visualization. | t2tv2.fasta |
| 3 | Coverage (2CN); *optional* | Expected diploid read coverage. | *unspecified* |
| 4 | Seed Overlap (bp); *optional* | Seed extension overlap threshold for support reads selection. | *unspecified* |
| 5 | Assembly Overlap (bp); *optional* | Extension overlap threshold for assembly. | *unspecified* |


### Output Files

| File | Description |
|---|---|
| *.metrics | Updated read statistics and coverage estimations for the final assembly. |
| *.ecDNA.fasta | The final reconstructed sequences of identified ecDNA structures. |
| *.refined.bam | Refined long-read alignments to the reconstructed ecDNA amplicons. |
| *.cycles.txt | A summary of the identified circular paths and segment coordinates. |
| *.ecLego3.run.log | The complete execution log, including timestamps and tool parameters. |

**TODO: write up the output files**


### Pipeline Limitations and Technical Refinements

- **Workflow Migration**: Transitioning the architecture to a Nextflow framework to facilitate dynamic resource allocation and optimize multi-core task parallelization.

- **Algorithmic Curation**: Implementation of rigorous automated filtering protocols to identify and exclude spurious circular genome assemblies.

- **Autonomous Stratification**: Development of fully automated pipelines for the characterization and distribution analysis of genetic subpopulations.

- **Modular Expansion**: Integration of comprehensive downstream modules for comparative analysis and functional annotation.


## Acknowledgements
* Our collaborators at Chang Gung Memorial Hospital, Linkou, Taiwan
* Albert H. Kim and the Tumor Bank of the Brain Tumor Center, St. Louis
* JGM GT team, esp. Chew Yee NGAN and Meihong LI
* Wei Lab; Aziz Taghbalout
* [Flye](https://github.com/fenderglass/Flye/) (https://github.com/fenderglass/Flye/).
* [KMC3](https://github.com/refresh-bio/KMC)(https://github.com/refresh-bio/KMC)
* [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)(https://github.com/fritzsedlazeck/Sniffles)


## Citations
"An Atlas of Extrachromosomal DNA Structures Illuminates Its Evolution and Biogenesis in Cancer" [DOI:10.64898/2025.12.24.696443](https://doi.org/10.64898/2025.12.24.696443)

