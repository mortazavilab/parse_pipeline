# Parse Biosciences Split-seq pre-processing workflow for HPC cluster
## Background: Combinatorial barcoding for snRNA-seq
Combinatorial or split-pool barcoding appends a unique **set** of barcodes to single cells or nuclei during multiple "rounds". These reactions are confined to each individual fixed and permeabilized cell, bypassing the need for physical barriers between cells such as microwells or nanodroplets.

This workflow is designed to pre-process mouse data from the [Parse Biosciences](https://www.parsebiosciences.com/) Evercode whole-transcriptome (“WT”) mini, 100k, and 1M barcoding kits. Versions 1, 2, and 3 of the barcoding chemistry are all supported. As of October 2024, genome and annotation is fixed at mm39 and vM32. This split-pool method may differ from others such as SHARE-seq due to using of both oligo-dT and random hexamer primers during the first round of reverse transcription-based barcoding. The goal in using both primers is to span the entire transcript for more full-length coverage and to avoid polyA-biased reads. After all 3 round of barcoding, aliquots of cells are lysed and libraries are prepared from barcoded cDNA. At the end of short-read library prep, a final Illumina barcode is added to identify these lysed aliquots, or "subpools", and libraries are sequenced and demultiplexed as subpools containing a mix of all the samples on the plate.

## Workflow overview

The main deviation from the official company pipeline is the mapping/quantification software and custom handling of sample merging to produce Python [anndata](https://anndata.readthedocs.io/en/latest/) objects grouped by tissue for downstream processing with [scanpy](https://scanpy.readthedocs.io/en/stable/). This workflow is based on the [kb-python](https://www.kallistobus.tools/kb_usage/kb_count/) (kallisto/bustools) suite, the official snRNA-seq tools for IGVF, and also includes optional code in place for genetic demultiplexing of individuals from sample barcoding wells loaded with 2 individuals with distinct genotypes.

<img align="left" src="https://github.com/mortazavilab/parse_pipeline/blob/main/uci_parse_workflow.png" width="422" height="451">

This workflow does the following:

1. **Pseudoalignment** with [kb count](https://github.com/pachterlab/kb_python): Associate reads with their cells/nuclei of origin, collapse duplicated reads according to unique molecular identifiers (UMIs), collapse reads from oligo-dT and random hexamer primers in the same well (`counts_unfiltered_modified`), and generate subpool-level cell x gene matrices.
2. **Genetic demultiplexing** with [klue](https://github.com/Yenaled/klue): Optionally, quantify reads associated with distinct genotypes for each cell & merge with cellbender adata. Assign genotype with higher count number between 2 expected genotypes, `tie` otherwise.
3. **Background removal** with [cellbender](https://cellbender.readthedocs.io/en/latest/index.html): Remove ambient RNA / PCR chimeras from data and filter cells/nuclei.
4. **Add metadata and doublet scores** with [scanpy](https://scanpy.readthedocs.io/en/stable/): Merge [sample metadata](https://github.com/mortazavilab/parse_pipeline/blob/main/configs/sample_metadata.csv) with cellbender output by sample barcode-to-well mapping to create a cellbender-filtered adata complete with sample and subpool information. Perform doublet scoring within scanpy using [scrublet](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html). 
5. **Merge** samples across subpools by tissue for downstream analysis.

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

## Setup instructions
Note: Written for and tested on UCI's HPC. Follow these steps if you are using this workflow for the first time.

### Clone this repository
Clone this branch.
```bash
git clone https://github.com/mortazavilab/parse_pipeline.git
```
### Set up an interactive tmux session on HPC
1. Remember your login node or choose your favorite out of i15, i16, i17. You can switch login nodes via ssh, e.g. `ssh login-i15`
2. Start tmux session via `tmux new -s mysession`, or whatever you want to name it. When you need to reconnect, type `tmux a -t mysession` ([tmux cheatsheet](https://tmuxcheatsheet.com/)) from the **same login node** as when you started the session. This is so if your internet goes out or you have to close your laptop, ongoing processes won't be terminated. 
3. Start interactive session: `srun -A SEYEDAM_LAB --cpus-per-task=1 --mem 32G --pty bash -i`. This is so you don't clog up the login node. We shouldn't need a lot since we are just installing packages and running snakemake, which will launch more computationally-intensive jobs for you.

### Create a conda environment called snakemake using python 3.9
1. `conda install -n base -c conda-forge mamba`
2. `mamba create -n snakemake -c conda-forge -c bioconda snakemake==7.32 python==3.9 scanpy==1.10.2 pandas==2.2.2 anndata==0.10.9 scikit-image==0.24.0`
3. `pip install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@04c2f5b460721fd55cf62a4cd23617b2555d69b8` for the latest CellBender version (3.2.0). (Yes unfortunately you need to install cellbender on this environment as well in order to use `anndata_from_h5` from `cellbender.remove_background.downstream` when making the adata from CellBender output.)

### Create a conda environment called cellbender using python 3.7
Tried multiple times to get CellBender working in the main snakemake environment with python 3.9 but to no avail; it's a [known version issue regarding saving the checkpoint file](https://github.com/broadinstitute/CellBender/issues/371). The checkpoint is very convenient in case the run fails due to lack of memory or job timeout, so I think the up-front effort in creating a second conda environment is worth it. Following the [installation instructions for CellBender](https://cellbender.readthedocs.io/en/latest/installation/index.html):

1. `mamba create -n cellbender python=3.7`
3. `mamba activate cellbender`
4. `mamba install -c anaconda pytables`
5. `pip install torch==1.13.1`
6. `pip install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@04c2f5b460721fd55cf62a4cd23617b2555d69b8` for the latest CellBender version (3.2.0).

### Install klue
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
2. `pwd` and copy the path, for example mine is `/dfs5/bio/erebboah/klue/build/src`
3. Then edit your `~/.bashrc` by adding this line to the bottom of the file: `export PATH="<your path>:$PATH"`. For example I added `export PATH="/dfs5/bio/erebboah/klue/build/src:$PATH`
4. `source ~/.bashrc`

Test installation by typing `klue` in the terminal, should see version and usage information.

## Create / edit input files
1. Make a fastq config file with columns that match the many examples in the configs folder, e.g. [igvf_004_config.tsv](https://github.com/mortazavilab/parse_pipeline/blob/cellbender_hpc/configs/igvf_004_config.tsv) and save in configs folder. Required columns are **fastq**, **fastq_r2**, **subpool**, **plate**, **lane**, **run**, **droplets_included**, **learning_rate**, and **expected_cells**. The last 3 parameters are specifically for CellBender [parameters of the same name](https://cellbender.readthedocs.io/en/latest/usage/index.html#recommended-best-practices). We recommend entering the same number of cells/nuclei as expected in the subpool (for example 67000 for 1M kits, 13000 for 100k kits) for expected cells, half the default learning rate (0.00005), and approximately 3 times the expected cells for the total droplets included (for example 200000 for 1M kits and 40000 for 100k kits). In the case of subpar CellBender results with these parameters, modify droplets included by including an additional 50,000 or subtracting 50,000 (1M kits).
2. Update [sample_metadata.csv](https://github.com/fairliereese/parse_pipeline/blob/main/configs/sample_metadata.csv) with relevant metadata for your experiment. The pipeline works on the assumption you have something in the following columns, in order to adhere to a common set of attributes especially when submitting data to the IGVF portal:
- Mouse_Tissue_ID
- Experiment
- bc1_well
- well_type
- SampleType
- Tissue
- Genotype
- Sex
- Age
- Protocol
- Chemistry
- Multiplexed_sample1
- Multiplexed_sample2
- Row
- Column
- DOB
- Age_days
- Body_weight_g
- Estrus_cycle
- Dissection_date
- Dissection_time
- ZT
- Dissector
- Tissue_weight_mg
- Notes
3. Make a copy of `snakemake/Snakefile.smk` with your experiment name to edit. Keep the copy in the snakemake folder. You can also just edit `Snakefile.smk` directly, but I like to keep a copy of the Snakefile used for each run for my own records. You only need to edit 4 lines maximum in the header region of the file:
- **config_tsv**: Path to the fastq config file which has the paths to your read1 and read2 input fastqs and the plate, lane, run, and sequencing platform.
- **sample_csv**: Path to the sample metadata. I typically update the [Google spreadsheet](https://docs.google.com/spreadsheets/d/13M6-Ry6oXgkx94BHZOGioYPI6F_hWDqjGgcaNu2JNYs/edit#gid=2080239896), download the tab, and upload it to my configs folder. Each row represents a well in the sample barcoding plate with metadata information for that sample, some of which is repeated across all the samples on the plate, such as experiment name, kit, and tissue.
- **kit**: either WT (48 wells), WT_mini (12 wells), or WT_mega (96 wells)
- **chemistry**: either v1, v2, or v3 

## Run pipeline
Skip steps 1-4 below if you were following the setup instructions and are already in an interactive tmux session.

1. Pay attention to your login node, or ssh to your favorite, e.g. `ssh login-i15`
2. Change directories to your pipeline directory, e.g. `cd /share/crsp/lab/seyedam/erebboah/parse_pipeline`. You MUST be in the `parse_pipeline` directory, not in a sub-directory like `parse_pipeline/configs`, `parse_pipeline/snakemake`, or it will not run.
3. Start tmux session, e.g. `tmux new -s mysession`
4. Start interactive session: `srun -A SEYEDAM_LAB --cpus-per-task=1 --mem 32G --pty bash -i`
5. Activate your snakemake environment: `conda activate snakemake` (you need to activate your snakemake conda environment again even if you were following setup instructions since sourcing the bashrc probably reset it to your base environment)
6. Check that snakemake is going to run the appropriate jobs (use the -n flag first). Make sure to change `Snakefile.smk` to the one you are actually using! For example, `Snakefile_igvf015.smk`
   
```bash
 snakemake \
  -s snakemake/Snakefile.smk \
  -j 50 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A {resources.account} --partition={resources.partition} --mem={resources.mem_gb}GB --gres={resources.gres} -c {resources.threads} --time=72:00:00" -n
 ```
7. Actually run pipeline (no `-n` flag):
```bash
 snakemake \
  -s snakemake/Snakefile.smk \
  -j 50 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A {resources.account} --partition={resources.partition} --mem={resources.mem_gb}GB --gres={resources.gres} -c {resources.threads} --time=72:00:00" 
 ```

### Basic troubleshooting
- FileNotFoundError/No such file or directory: Check your current directory (`pwd`). Make sure the 3 required input files exist and in the correct locations: fastq config e.g. `igvf_###_config.tsv` is in `parse_pipeline/configs`, `sample_metadata.csv` is in `parse_pipeline/configs`, and `Snakemake_###.smk` is in `parse_pipeline/snakemake`. Make sure the fastq config file is spelled correctly in your Snakemake smk file.
- AttributeError: Make sure the columns in `igvf_###_config.tsv` exactly match **fastq**, **fastq_r2**, **subpool**, **plate**, **lane**, **run**, **droplets_included**, **learning_rate**, and **expected_cells**.
- If your fastq or fastq_r2 paths are incorrect, the pipeline will run fine at first. Kallisto runs, symlinked files have the correct name while being wrong, but you get no cells and code downstream of kallisto will probably break.
- If you put the wrong kit/chemistry in the snakemake file, you will get NAN barcodes and the pipeline fails due to an assert() call in the make cellbender adata rule.

### Known issues / Wishlist
- Clean up extra files
- klue reference generation runs twice for the F1 plates, in other words it makes both NODJ.idx file and a B6NODF1J.idx file, which for now are identical. And the same for all the other 6 non-B6J genotypes. The good news is that it only runs twice one time...reference generation isn't repeated after the first pipeline run.
- Integrate report code
- Integrate experimentSpecs and analysisSpecs to replace manual metadata curation
- Support human genome
- Support custom humanized loci for Model AD
