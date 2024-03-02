m# Parse Biosciences Split-seq Pre-processing Workflow
## Background: Combinatorial barcoding for snRNA-seq
Combinatorial or split-pool barcoding appends a unique **set** of barcodes to single cells or nuclei during multiple "rounds". These reactions are confined to each individual fixed and permeabilized cell, bypassing the need for physical barriers between cells such as microwells or nanodroplets.

This workflow is designed to pre-process data from the [Parse Biosciences](https://www.parsebiosciences.com/) Evercode whole-transcriptome (“WT”) barcoding kits. This split-pool method may differ from others such as [SHARE-seq](https://www.sciencedirect.com/science/article/pii/S0092867420312538?via%3Dihub) due to using of both oligo-dT and random hexamer primers during the first round of reverse transcription-based barcoding. The goal in using both primers is to span the entire transcript for more full-length coverage and to avoid polyA-biased reads. After all 3 round of barcoding, aliquots of cells are lysed and libraries are prepared from barcoded cDNA. At the end of short-read library prep, a final Illumina barcode is added to identify these lysed aliquots, or "subpools", and libraries are sequenced and demultiplexed as subpools containing a mix of all the samples on the plate.

## Workflow overview

The main deviation from the official company pipeline is the mapping/quantification software and custom handling of sample merging to produce Python [anndata](https://anndata.readthedocs.io/en/latest/) objects grouped by tissue for downstream processing with [scanpy](https://scanpy.readthedocs.io/en/stable/). This workflow is based on the [kb-python](https://www.kallistobus.tools/kb_usage/kb_count/) (kallisto/bustools) suite, the official snRNA-seq tools for IGVF, and also includes optional code in place for genetic demultiplexing of individuals from sample barcoding wells loaded with 2 individuals with distinct genotypes.

This workflow does the following:

1. **[kb-count](https://github.com/pachterlab/kb_python)**: Associate reads with their cells/nuclei of origin
2. **kb-count**: Collapse duplicated reads according to unique molecular identifiers (UMIs)
3. **kb-count**: Collapse reads from oligo-dT and random hexamer primers in the same well (`counts_unfiltered_modified`)
4. **kb-count**: Generate cell x gene adatas for each subpool
5. **[klue](https://github.com/Yenaled/klue)**: Quantify reads associated with distinct genotypes for each cell & merge with kallisto adata
6. Merge custom metadata by sample barcode-to-well mapping
7. Assign genotype with higher count number between 2 expected genotypes, `tie` otherwise 
8. Perform basic filtering for easier merging (adjust in config.yml, default 200 UMI)
9. **[scrublet](https://github.com/swolock/scrublet)**: Generate doublet scores 
10. Merge samples across subpools by tissue of origin

## Input data
Subpool FASTQ Files:
- Subpool_1_R1.fastq.gz (contains cDNA)
- Subpool_1_R2.fastq.gz (contains the cell barcode and UMI)

Fastq files split by multiple lanes are ok, e.g. 
- Subpool_2_S1_L001_R1_001.fastq.gz
- Subpool_2_S1_L001_R2_001.fastq.gz
- Subpool_2_S1_L003_R1_001.fastq.gz
- Subpool_2_S1_L003_R2_001.fastq.gz
etc.

Read name is not used directly in the pipeline, can be formatted however. Just need to specify the relevant information in the fastq config file.

# Setup instructions
Note: Written for and tested on UCI's HPC cluster with slurm job scheduler. Follow these steps if you are using this workflow for the first time.

## Clone this repository
Choose a location on HPC with plenty of space to clone this repo. The kallisto output in particular is very large for a full 1M cell (Mega kit) experiment (~**350GB** with all unfiltered bus files).
```bash
git clone https://github.com/fairliereese/parse_pipeline.git
```
## Set up an interactive tmux session on HPC
1. Remember your login node or choose your favorite out of i15, i16, i17. You can switch login nodes via ssh, e.g. `ssh login-i15`
2. Start tmux session via `tmux new -s mysession`, or whatever you want to name it. When you need to reconnect, type `tmux a -t mysession` ([tmux cheatsheet](https://tmuxcheatsheet.com/)) from the **same login node** as when you started the session. This is so if your internet goes out or you have to close your laptop, ongoing processes won't be terminated. 
3. Start interactive session: `srun -A SEYEDAM_LAB --cpus-per-task=1 --mem 32G --pty bash -i`. This is so you don't clog up the login node. If you ever see the error 'Killed', you are probably on the login node, or you need to increase the requested memory for the interactive session. We shouldn't need a lot since we are just installing packages and running snakemake, which will launch more computationally-intensive jobs for you.

## Create a conda environment called snakemake
Required packages: `snakemake`, `pandas`, `numpy`, `anndata`, `scanpy`, `scrublet`, `kb-python`, and if you have genetically multiplexed samples, `klue`.
1. `conda install -n base -c conda-forge mamba`
2. `conda create -c conda-forge -c bioconda -n snakemake snakemake==7.32 python==3.9 pandas`
3. `conda activate snakemake`
4. Install required python packages with pip,
`pip install kb-python scrublet`

Klue installation instructions:
1. `git clone https://github.com/Yenaled/klue`
2. `cd klue`
3. `mkdir build`
4. `cd build`
5. `cmake ..`
6. `make`
7. `make install`

If step 7 does not work, follow these steps:
1. `cd src` (within the build folder)
2. `pwd` and copy the path, for example mine is `/share/crsp/lab/seyedam/erebboah/parse_pipeline/klue/build/src`
3. Then edit your `~/.bashrc` by adding this line to the bottom of the file: `export PATH="<your path>:$PATH"`. For example I added `export PATH="/share/crsp/lab/seyedam/erebboah/parse_pipeline/klue/build/src:$PATH"`
4. `source ~/.bashrc`

Test installation by typing `klue` in the terminal, should see version and usage information.

# Create / update supplemental files
1. Fastq config file e.g. [igvf_003_config.tsv](https://github.com/fairliereese/parse_pipeline/blob/main/configs/igvf_003_config.tsv) and save in configs folder
2. Update [sample metadata file](https://github.com/fairliereese/parse_pipeline/blob/main/configs/sample_metadata.csv) with relevant metadata for your experiment. The minimum required metadata columns for the pipeline to run properly are **Mouse_Tissue_ID**, **Experiment**, **bc1_well**, **well_type**, **Tissue**, and **Genotype**. If you have genetically multiplexed wells, it's also very convienent for downstream analysis to also have **Multiplexed_sample1** and **Multiplexed_sample2**.
3. Update snakemake file with name of your fastq config and check to make sure kit and chemistry are correct. - TODO make example

# Run pipeline
Skip steps 1-4 if you were following the setup instructions and are already in an interactive tmux session. 

1. Pay attention to your login node, or ssh to your favorite, e.g. `ssh login-i15`
2. Change directories to the your pipeline directory, e.g. `cd /share/crsp/lab/seyedam/erebboah/parse_pipeline`
3. Start tmux session, e.g. `tmux new -s mysession`
4. Start interactive session: `srun -A SEYEDAM_LAB --cpus-per-task=1 --mem 32G --pty bash -i`
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


