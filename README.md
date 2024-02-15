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
Set up fastq specification file and move to here `/share/crsp/lab/seyedam/share/igvf_pipeline/configs`

## Run pipeline
1. Pay attention to your login node — or choose your favorite out of i15, i16, i17. `ssh login-i15`
2. go here: /share/crsp/lab/seyedam/share/igvf_pipeline and
3. start tmux session via tmux new -s mysession  If you need to reconnect — tmux a -t mysession (https://tmuxcheatsheet.com/).
4. start interactive session srun -A SEYEDAM_LAB --cpus-per-task=1 --time=168:00:00 --mem 8GB --pty bash -i
5. activate your snakemake environment: conda activate snakemake
6. Check that snakemake is going to run the appropriate jobs (use the -n flag first)
```bash
 snakemake \
  -s snakemake/Snakefile_016 \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --time=72:00:00" -n
 ```
7. Actually run pipeline
```bash
snakemake   -s snakemake/Snakefile_016   -j 100   --latency-wait 120   --use-conda   --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --time=72:00:00"
 ```
