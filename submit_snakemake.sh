#!/bin/bash

#SBATCH -J wgbs_snakemake         # Job name
#SBATCH -p batch                  # Partition
#SBATCH --qos long_jobs           # long jobs
#SBATCH -N 1                      # Number of nodes
#SBATCH --cpus-per-task=4         # cpu per task
#SBATCH --mem=8G                  # memory
#SBATCH -t 4-00                   # Maximum runtime
#SBATCH -o logs/snakemake_%j.out  # Standard output log
#SBATCH -e logs/snakemake_%j.err  # Standard error log

source activate /home/groups/hoolock2/u0/bd/miniconda3/envs/snake

snakemake --use-conda --jobs 50 --latency-wait 60 --keep-going --rerun-incomplete --cluster-config cluster.json --cluster "sbatch --qos {cluster.qos} -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem {cluster.mem} -t {cluster.time} -o {cluster.stdout} -e {cluster.stderr} -J {rule}.{wildcards}"
