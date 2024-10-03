# Parse Biosciences Split-seq pre-processing workflow for HPC cluster
## Background: Combinatorial barcoding for snRNA-seq
Combinatorial or split-pool barcoding appends a unique **set** of barcodes to single cells or nuclei during multiple "rounds". These reactions are confined to each individual fixed and permeabilized cell, bypassing the need for physical barriers between cells such as microwells or nanodroplets.

This workflow is designed to pre-process mouse data from the [Parse Biosciences](https://www.parsebiosciences.com/) Evercode whole-transcriptome (“WT”) mini, 100k, and 1M barcoding kits. Versions 1, 2, and 3 of the barcoding chemistry are all supported. As of October 2024, genome and annotation is fixed at mm39 and vM32. This split-pool method may differ from others such as [SHARE-seq](https://www.sciencedirect.com/science/article/pii/S0092867420312538?via%3Dihub) due to using of both oligo-dT and random hexamer primers during the first round of reverse transcription-based barcoding. The goal in using both primers is to span the entire transcript for more full-length coverage and to avoid polyA-biased reads. After all 3 round of barcoding, aliquots of cells are lysed and libraries are prepared from barcoded cDNA. At the end of short-read library prep, a final Illumina barcode is added to identify these lysed aliquots, or "subpools", and libraries are sequenced and demultiplexed as subpools containing a mix of all the samples on the plate.

## Workflow overview

The main deviation from the official company pipeline is the mapping/quantification software and custom handling of sample merging to produce Python [anndata](https://anndata.readthedocs.io/en/latest/) objects grouped by tissue for downstream processing with [scanpy](https://scanpy.readthedocs.io/en/stable/). This workflow is based on the [kb-python](https://www.kallistobus.tools/kb_usage/kb_count/) (kallisto/bustools) suite, the official snRNA-seq tools for IGVF, and also includes optional code in place for genetic demultiplexing of individuals from sample barcoding wells loaded with 2 individuals with distinct genotypes.

This workflow does the following:

1. **[kb-count](https://github.com/pachterlab/kb_python)**: Associate reads with their cells/nuclei of origin
2. **kb-count**: Collapse duplicated reads according to unique molecular identifiers (UMIs)
3. **kb-count**: Collapse reads from oligo-dT and random hexamer primers in the same well (`counts_unfiltered_modified`)
4. **kb-count**: Generate an unfiltered cell x gene adata for each subpool
5. (HOPEFULLY) **[cellbender](https://cellbender.readthedocs.io/en/latest/index.html)**: Remove ambient RNA / PCR chimeras from data and filter cells/nuclei
6. **scanpy**: Merge [sample metadata](https://github.com/mortazavilab/parse_pipeline/blob/main/configs/sample_metadata.csv) with cellbender output by sample barcode-to-well mapping to create a cellbender-filtered adata complete with sample and subpool information
7. **[klue](https://github.com/Yenaled/klue)**: Quantify reads associated with distinct genotypes for each cell & merge with cellbender adata. Assign genotype with higher count number between 2 expected genotypes, `tie` otherwise 
8. **[scrublet](https://github.com/swolock/scrublet)**: Generate doublet scores 
9. Merge samples across subpools by genotype & tissue of origin

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
Note: Written for and tested on the Mortazavi lab's Watson server. Follow these steps if you are using this workflow for the first time.

### Clone this repository
Clone this repo. The kallisto output in particular is very large for a full 1M cell (Mega kit) experiment (~**350GB** with all unfiltered bus files).
```bash
git clone -b cellbender https://github.com/mortazavilab/parse_pipeline.git
```
### Set up a screen session on Watson
Start screen session via `screen -S mysession`, or whatever you want to name it. When you need to reconnect, type `screen -r mysession` ([screen cheatsheet](https://kapeli.com/cheat_sheets/screen.docset/Contents/Resources/Documents/index)). This is so if your internet goes out or you have to close your laptop, ongoing processes won't be terminated. 

### Create a conda environment called snakemake
Required packages: `snakemake`, `pandas`, `numpy`, `anndata`, `scanpy`, `scrublet`, `kb-python`, and if you have genetically multiplexed samples, `klue`.
1. `conda install -n base -c conda-forge mamba`
2. `mamba create -c conda-forge -c bioconda -n snakemake snakemake==7.32 python==3.9 pandas`
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
2. `pwd` and copy the path, for example mine is `/home/erebboah/klue/build/src`
3. Then edit your `~/.bashrc` by adding this line to the bottom of the file: `export PATH="<your path>:$PATH"`. For example I added `export PATH="/home/erebboah/klue/build/src:$PATH"`
4. `source ~/.bashrc`

Test installation by typing `klue` in the terminal, should see version and usage information.

## Create / edit input files
1. Make a fastq config file with columns that match the many examples in the configs folder, e.g. [igvf_015_config.tsv](https://github.com/fairliereese/parse_pipeline/blob/main/configs/igvf_015_config.tsv) and save in configs folder. Required columns are **fastq**, **fastq_r2**, **subpool**, **plate**, **lane**, and **run**.
2. Update [sample_metadata.csv](https://github.com/fairliereese/parse_pipeline/blob/main/configs/sample_metadata.csv) with relevant metadata for your experiment. The minimum required metadata columns for the pipeline to run properly are **Mouse_Tissue_ID**, **Experiment**, **bc1_well**, **well_type**, **Tissue**, and **Genotype**. If you have genetically multiplexed wells, it's also very convienent for downstream analysis to also have **Multiplexed_sample1** and **Multiplexed_sample2**.
3. Make a copy of `snakemake/Snakefile.smk` with your experiment name to edit. Keep the copy in the snakemake folder. You can also just edit `Snakefile.smk` directly, but I like to keep a copy of the Snakefile used for each run for my own records. You only need to edit 5 lines maximum in the header region of the file:
- **config_tsv**: Path to the fastq config file which has the paths to your read1 and read2 input fastqs and the plate, lane, run, and sequencing platform.
- **sample_csv**: Path to the sample metadata. I typically update the [Google spreadsheet](https://docs.google.com/spreadsheets/d/13M6-Ry6oXgkx94BHZOGioYPI6F_hWDqjGgcaNu2JNYs/edit#gid=2080239896), download the tab, and upload it to my configs folder. Each row represents a well in the sample barcoding plate with metadata information for that sample, some of which is repeated across all the samples on the plate, such as experiment name, kit, and tissue.
- **kit**: either WT (48 wells), WT_mini (12 wells), or WT_mega (96 wells)
- **chemistry**: all IGVF and ModelAD experiments are v2 so far, v3 coming soon
- **first_min_counts**: minimum # of UMIs per cell
 

## Run pipeline
Skip steps 1-4 if you were following the setup instructions and are already in an interactive tmux session.

1. Login to Watson and start a screen session, e.g. `screen -S pipe`
2. Change directories to your pipeline directory, e.g. `cd /home/erebboah/parse_pipeline`. You MUST be in the `parse_pipeline` directory, not in a sub-directory like `parse_pipeline/configs`, `parse_pipeline/snakemake`, or it will not run.
3. Activate your snakemake environment: `conda activate snakemake` (you need to activate your snakemake conda environment again even if you were following setup instructions since sourcing the bashrc probably reset it to your base environment)
6. Check that snakemake is going to run the appropriate jobs (use the -n flag first). Make sure to change `Snakefile.smk` to the one you are actually using! For example, `Snakefile_igvf015.smk`
   
```bash
snakemake  -s snakemake/Snakefile.smk --latency-wait 120  --use-conda --cores 1 -n
 ```
7. Actually run pipeline
```bash
snakemake  -s snakemake/Snakefile.smk --latency-wait 120  --use-conda --cores 1
 ```

### Basic troubleshooting
- FileNotFoundError/No such file or directory: Check your current directory (`pwd`). Make sure the 3 required input files exist and in the correct locations: fastq config e.g. `igvf_###_config.tsv` is in `parse_pipeline/configs`, `sample_metadata.csv` is in `parse_pipeline/configs`, and `Snakemake_###.smk` is in `parse_pipeline/snakemake`. Make sure the fastq config file is spelled correctly in your Snakemake smk file.
- AttributeError: Make sure the columns in `igvf_###_config.tsv` exactly match **fastq**, **fastq_r2**, **subpool**, **plate**, **lane**, and **run**.

### Known issues / Wishlist
- CELLBENDER PLZ
- Clean up extra files
- Scrublet woes...
- klue reference generation runs twice for the F1 plates, in other words it makes both NODJ.idx file and a B6NODF1J.idx file, which for now are identical. And the same for all the other 6 non-B6J genotypes. The good news is that it only runs twice one time...reference generation isn't repeated after the first pipeline run.
- Integrate Ryan's report code to make beautiful knee plots and well heatmaps
- Integrate experimentSpecs and analysisSpecs to replace manual metadata curation
- Support custom humanized loci for Model AD
