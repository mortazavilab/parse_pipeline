```
conda activate base_clone
snakemake \
  -s snakemake/Snakefile \
  -j 100 \
  --latency-wait 120 \
  -n
  ```

 ```bash
 snakemake \
  -s snakemake/Snakefile \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
 ```
