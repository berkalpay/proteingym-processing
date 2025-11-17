The script [extract_metadata.py](scripts/extract_metadata.py) is used to first extract the collected metadata file into separate single-entry CSVs, allowing Snakemake to detect whether there have been any updates to the metadata.
The processing pipeline is expressed in [Snakefile](Snakefile): an initial alignment is made and the bitscore is revised until certain criteria are met.
See [example.sbatch](example.sbatch) for an example of how to use snakemake to make a finalized alignment.
Additional steps such as score processing are not yet fully implemented.
