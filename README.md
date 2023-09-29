```
conda activate base_clone
snakemake \
  -s snakemake/Snakefile \
  -j 100 \
  --latency-wait 120 \
  -n
  ```
