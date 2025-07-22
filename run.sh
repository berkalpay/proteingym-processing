snakemake -n -p --use-conda \
  -j 100 \
  --cluster-config cluster-config.json \
  --cluster "sbatch --mem={cluster.mem} -t {cluster.time} -c {cluster.cores} -p {cluster.partition} -G {cluster.gpu} --output=logs-jobs/slurm-%j.out" \
  --keep-going \
  --latency-wait 60 \
  --restart-times 0