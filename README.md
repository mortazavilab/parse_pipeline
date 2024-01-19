```bash
 snakemake \
  -s snakemake/Snakefile \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --time=72:00:00" -n
 ```


 ```bash
conda activate snakemake_vis
snakemake -s snakemake/Snakefile --forceall --dag | dot -Tpdf > dag.pdf
snakemake -s snakemake/Snakefile --forceall --rulegraph | dot -Tpdf > dag.pdf
```
