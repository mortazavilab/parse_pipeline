#!/bin/sh
#SBATCH --job-name=scpp   ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH --nodes=1 ## (-N) number of nodes to use
#SBATCH -p highmem
#SBATCH --cpus-per-task=16         ## number of cores the job needsi
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err
#SBATCH --mem=256G
#SBATCH --array=1-9

source /data/homezvol1/weberrl/miniconda3/bin/activate singlecell

tissue_files=('Adrenalpreprocessed.h5ad' 'Gastrocnemiuspreprocessed.h5ad' 'GonadsMalepreprocessed.h5ad' 'HypothalamusPituitarypreprocessed.h5ad' 'Liverpreprocessed.h5ad' 'CortexHippocampuspreprocessed.h5ad' 'GonadsFemalepreprocessed.h5ad' 'Heartpreprocessed.h5ad' 'Kidneypreprocessed.h5ad')

# Get the index corresponding to the array job task ID
index=$((SLURM_ARRAY_TASK_ID - 1))

# Extract the filename using the index
filename=${tissue_files[$index]}

# Run the Python script with the selected filename
python scanpy_preprocessing.py -f $filename
