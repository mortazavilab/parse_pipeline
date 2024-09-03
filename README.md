# Parse Biosciences Split-seq pre-processing workflow for local machine with GPU

## Workflow overview

Starting from the raw kallisto output and if applicable the klue genotype counts, this workflow does the following:

1. [cellbender](https://cellbender.readthedocs.io/en/latest/index.html)**: Remove ambient RNA / PCR chimeras from data and filter cells/nuclei
2. **scanpy**: Merge [sample metadata](https://github.com/mortazavilab/parse_pipeline/blob/main/configs/sample_metadata.csv) with cellbender output by sample barcode-to-well mapping to create a cellbender-filtered adata complete with sample and subpool information
3. **[scrublet]([https://github.com/swolock/scrublet](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html))**: Calculate doublet scores using scrublet called from within scanpy.
4. Merge samples across subpools for a plate-level adata.

## Input data
1. Raw kallisto adata
2. Klue genotype count tsv file

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
Skip steps 1-4 if you were following the setup instructions and are already in an interactive screen session.

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
- klue reference generation runs twice for the F1 plates, in other words it makes both NODJ.idx file and a B6NODF1J.idx file, which for now are identical. And the same for all the other 6 non-B6J genotypes. The good news is that it only runs twice one time...reference generation isn't repeated after the first pipeline run.
- Integrate report code
- Integrate experimentSpecs and analysisSpecs to replace manual metadata curation
- Support custom humanized loci for Model AD
