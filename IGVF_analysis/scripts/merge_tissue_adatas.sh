#!/bin/sh
#SBATCH --job-name=merge   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH --nodes=1 ## (-N) number of nodes to use
#SBATCH -p highmem
#SBATCH --cpus-per-task=16         ## number of cores the job needsi
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err
#SBATCH --mem=256G

source /data/homezvol1/weberrl/miniconda3/bin/activate singlecell


python merge_tissue_adatas.py
