rule extract_metadata:
    input:
        ancient("data/DMS_substitutions_internal.csv")
        "{dms}"
    output:
        "data/metadata/{scan}.csv"
    shell:
        "python scripts/extract_metadata.py {input} > {output}"  # only make if new?
    conda:
        "envs/base.yaml"


rule align:
    input:
        "data/metadata/{scan}.csv"
    output:
        "data/alignments/{scan}.a2m"
    shell:
        "python scripts/align.py {input} > {output}"
    conda:
        "envs/evcouplings.yaml"


rule process_scores:
    input:
        "data/metadata/{scan}.csv"
        "data/raw_scores/{scan}.csv"
    output:
        "data/scores/{scan}.csv"
        "data/binarization_cutoffs/{scan}.txt"
    shell:
        "python scripts/process_scores.py {input} > {output}"
    conda:
        "envs/base.yaml"
