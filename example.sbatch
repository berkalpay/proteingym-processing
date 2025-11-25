#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -t 9-00:00
#SBATCH -p long
#SBATCH --mem=3GB
#SBATCH --output=logs/slurm-%j.out


python3 scripts/extract_metadata.py data/DMS_substitutions_internal.csv data/metadata/

# NOTE: make sure to activate an environment that has snakemake
# Tells snakemake to generate the final alignment for this protein
snakemake data/alignments/DYR_ECOLI_Thompson_2019.a2m -p --use-conda --profile slurm --retries 15
