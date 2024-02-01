#!/bin/bash
#SBATCH --job-name=degs    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=degs-%J.out ## output log file
#SBATCH --error=degs-%J.err ## error log file
#SBATCH --mem=128G
#SBATCH --time=3-00:00

source ~/miniconda3/bin/activate hpc3sc

mkdir degs/Gastrocnemius

python3 pseudobulk_pydeseq.py --tissue Gastrocnemius
