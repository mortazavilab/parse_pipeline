#!/bin/bash
#SBATCH --job-name=kallisto    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=8         ## number of cores the job needs
#SBATCH --output=kallisto-%J.out ## output log file
#SBATCH --error=kallisto-%J.err ## error log file
#SBATCH --mem=64G
#SBATCH --array=1-15



source ~/miniconda3/bin/activate kb
igvf='igvf_010'
#igvf=$1
mkdir ../${igvf}
mkdir /dfs7/samlab/seyedam/IGVF/${igvf}/temp/


mkdir ../${igvf}/aj_pwk



all_sublibraries=('2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16')
samples=('S1' 'S2' 'S3' 'S4' 'S5' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15')

sublibraries=${all_sublibraries[$SLURM_ARRAY_TASK_ID - 1]}
current_sample=${samples[$SLURM_ARRAY_TASK_ID - 1]}

nova1_path='/dfs7/samlab/seyedam/IGVF/igvf_010/nova1/'
nova2_path='/dfs7/samlab/seyedam/IGVF/igvf_010/nova2/'

ref_path='/share/crsp/lab/seyedam/share/igvf_splitseq/ref'

t2g='/share/crsp/lab/seyedam/share/igvf_splitseq/klue_demultiplexing/idx_t2g_files/distinguish.Mus_musculus.Mus_musculus_aj_Mus_musculus_pwkphj.everything.t2g'
idx='/share/crsp/lab/seyedam/share/igvf_splitseq/klue_demultiplexing/idx_t2g_files/distinguish.Mus_musculus.Mus_musculus_aj_Mus_musculus_pwkphj.everything.idx'




kb count --h5ad --gene-names --sum=nucleus --strand=forward -r ${ref_path}/r1_RT_replace.txt -w ${ref_path}/r1r2r3.txt --workflow=standard -g ${t2g} -x SPLIT-SEQ -i ${idx} -t 24 -o ../${igvf}/aj_pwk/${igvf}_sub${sublibraries} ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L001_R1_001.fastq.gz ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L001_R2_001.fastq.gz ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L002_R1_001.fastq.gz ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L002_R2_001.fastq.gz ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L003_R1_001.fastq.gz ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L003_R2_001.fastq.gz ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L004_R1_001.fastq.gz ${nova1_path}Sublibrary_${sublibraries}_${current_sample}_L004_R2_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L001_R1_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L001_R2_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L002_R1_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L002_R2_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L003_R1_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L003_R2_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L004_R1_001.fastq.gz ${nova2_path}Sublibrary_${sublibraries}_${current_sample}_L004_R2_001.fastq.gz


