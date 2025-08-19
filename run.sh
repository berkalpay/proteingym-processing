python3 scripts/extract_metadata.py data/DMS_substitutions_internal.csv data/metadata/

snakemake -n -p --use-conda --profile slurm
