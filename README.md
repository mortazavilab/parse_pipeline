# Parse Biosciences Split-seq
## Overview
Split-pool barcoding appends a unique **set** of barcodes to single cells or nuclei during multiple "rounds". These reactions are confined to each individual fixed and permeabilized cell, bypassing the need for physical barriers between cells such as microwells or nanodroplets.

This workflow is designed to pre-process data from the [Parse Biosciences](https://www.parsebiosciences.com/) Evercode whole-transcriptome (“WT”) barcoding kits. This split-pool method may differ from others such as SHARE-seq due to using of both oligo-dT and random hexamer primers during the first round of reverse transcription-based barcoding. The goal in using both primers is to span the entire transcript for more full-length coverage and to avoid polyA-biased reads. After all 3 round of barcoding, a final “subpool” generation step adds the final Illumina barcode and libraries are sequenced and demultiplexed as subpools containing a mix of all the samples on the plate.

The main deviation from the official company software is the mapping/quantification software and custom handling of sample merging to produce Python [anndata](https://anndata.readthedocs.io/en/latest/) objects grouped by tissue. This workflow is based on the [kb-python](https://www.kallistobus.tools/kb_usage/kb_count/) (kallisto/bustools) suite of tools, the official snRNA-seq software for IGVF, and also includes optional code in place for genetic demultiplexing using [klue](https://github.com/Yenaled/klue) of individuals from barcoding wells loaded with 2 individuals with distinct genotypes.

This workflow does the following:

1. **[kb-count](https://github.com/pachterlab/kb_python)**: Associate reads with their cells/nuclei of origin
2. **[kb-count](https://github.com/pachterlab/kb_python)**: Collapse duplicated reads according to unique molecular identifiers (UMIs)
3. **[kb-count](https://github.com/pachterlab/kb_python)**: Collapse reads from oligo-dT and random hexamer primers in the same well (`counts_unfiltered_modified`)
4. **[kb-count](https://github.com/pachterlab/kb_python)**: Generate cell x gene adatas for each subpool
5. **[klue](https://github.com/Yenaled/klue)**: Quantify reads associated with distinct genotypes for each cell & merge with kallisto adata
6. Merge custom metadata by sample barcode-to-well mapping
7. Assign genotype with higher count number between 2 expected genotypes, `tie` otherwise 
8. Perform basic filtering for easier merging (adjust in config.yml, default 200 UMI)
9. **[scrublet](https://github.com/swolock/scrublet)**: Generate doublet scores 
10. Merge samples across subpools by tissue of origin

## Input data
Subpool FASTQ Files:
Subpool_1_R1.fastq.gz (contains cDNA)
Subpool_1_R2.fastq.gz (contains the cell barcode and UMI)

# Instructions
Note: Written for and tested on UCI's HPC cluster. 

## Create a conda environment called snakemake
Required packages: `snakemake`, `pandas`, `numpy`, `anndata`, `scanpy`, `scrublet`, `kb-python`, and if you have genetically multiplexed samples, `klue`.
1. `conda install -n base -c conda-forge mamba`
2. `conda create -c conda-forge -c bioconda -n snakemake snakemake==7.32 python==3.9 pandas` 
3. Install required python packages with pip, e.g.
```bash
pip install kb-python
pip install scrublet
 ```

## Create required files
Set up fastq specification file e.g. [igvf_003_config.tsv](https://github.com/fairliereese/parse_pipeline/blob/main/configs/igvf_003_config.tsv) and move to here `/share/crsp/lab/seyedam/share/igvf_pipeline/configs`

## Run pipeline
1. Pay attention to your login node — or choose your favorite out of i15, i16, i17. `ssh login-i15`
2. Change directories to the main pipeline directory: `cd /share/crsp/lab/seyedam/share/igvf_pipeline`
3. Start tmux session via `tmux new -s mysession`  If you need to reconnect — `tmux a -t mysession` ([tmux cheatsheet](https://tmuxcheatsheet.com/)).
4. Start interactive session `srun -A SEYEDAM_LAB --cpus-per-task=1 --time=168:00:00 --mem 8GB --pty bash -i`
5. Activate your snakemake environment: `conda activate snakemake`
6. Check that snakemake is going to run the appropriate jobs (use the -n flag first)
   
```bash
 snakemake \
  -s snakemake/Snakefile.smk \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --time=72:00:00" -n
 ```
7. Actually run pipeline
```bash
snakemake \
-s snakemake/Snakefile.smk \
-j 100 \
--latency-wait 120 \
--use-conda \
--cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --time=72:00:00"
 ```
