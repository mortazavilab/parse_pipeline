#!/bin/sh
#SBATCH --job-name=scpp   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH --nodes=1 ## (-N) number of nodes to use
#SBATCH -p highmem
#SBATCH --cpus-per-task=1         ## number of cores the job needsi
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mem=256G
#SBATCH --array=1

source ~/miniconda3/bin/activate scanpy_env

tissue_files=('PBMC.h5ad')

# Get the index corresponding to the array job task ID
index=$((SLURM_ARRAY_TASK_ID - 1))

# Extract the filename using the index
filename=${tissue_files[$index]}

# Run the Python script with the selected filename
python scanpy_preprocessing_liz.py -f $filename
