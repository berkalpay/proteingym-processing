from glob import glob
from snakemake.io import glob_wildcards

(scans,) = glob_wildcards("data/metadata/{scan}.csv")


rule all:
    input:
        expand("data/alignments/{scan}.a2m", scan=scans),


def alignment_decision(summary_filepath: str) -> float | None:
    import os
    import pandas as pd

    bitscore = float(
        summary_filepath.split(os.sep)[-1].split("_alignment_statistics.csv")[0]
    )
    STEP = 0.05

    alignment_stats = pd.read_csv(summary_filepath)
    perc_cov = alignment_stats.at[0, "perc_cov"]
    num_seqs = alignment_stats.at[0, "num_seqs"]
    neff_over_l = alignment_stats.at[0, "N_eff"] / alignment_stats.at[0, "seqlen"]

    if (perc_cov < 0.7 or num_seqs < 100 or neff_over_l < 1) and bitscore - STEP > 0:
        return round(bitscore - STEP, 2)
    if neff_over_l > 100 and bitscore + STEP < 2:
        new_bitscore = round(bitscore + STEP, 2)
        return round(bitscore + 2 * STEP, 2) if new_bitscore == 1 else new_bitscore
    return None


def stats(wildcards):
    import os
    import pandas as pd

    alignments_dir = f"data/EVcouplings/{wildcards.scan}"
    stats_files = []
    if os.path.isdir(alignments_dir):  # Read stats of existing alignments
        for bitscore in os.listdir(alignments_dir):
            bitscore_dir = os.path.join(alignments_dir, bitscore, "align")
            if not os.path.isdir(bitscore_dir):
                continue
            for f in os.listdir(bitscore_dir):
                if f.endswith("_alignment_statistics.csv"):
                    stats_files.append(os.path.join(bitscore_dir, f))

    if not stats_files:  # Make initial alignment
        metadata = pd.read_csv(f"data/metadata/{wildcards.scan}.csv")
        initial_bitscore = 0.2 if metadata.at[0, "taxon"] == "Virus" else 0.7
        a2m = checkpoints.align.get(
            scan=wildcards.scan, bitscore=initial_bitscore
        ).output[0]
        return a2m.replace(".a2m", "_alignment_statistics.csv")

    return max(stats_files, key=lambda p: os.stat(p).st_mtime)


rule align:
    input:
        stats,
    output:
        "data/alignments/{scan}.a2m",
    conda:
        "/n/groups/marks/projects/marks_lab_and_oatml/ProteinGym2/EVCouplings/envs/pg2_evc"
    run:
        import os
        import shutil

        stats_file = str(input[0])
        new_bitscore = alignment_decision(stats_file)

        if new_bitscore is None:  # Finalize alignment
            shutil.copy(
                stats_file.replace("_alignment_statistics.csv", ".a2m"), output[0]
            )
        else:
            os.system(
                f"python scripts/align.py data/metadata/{wildcards.scan}.csv {new_bitscore}"
            )
            raise ValueError(
                f"Triggered new alignment at bitscore={new_bitscore}, "
                f"producing {a2m}; re-run needed to continue."
            )  # Make Snakemake rerun this rule


rule process_scores:
    input:
        "data/metadata/{scan}.csv" "data/raw_scores/{scan}.csv",
    output:
        "data/scores/{scan}.csv" "data/binarization_cutoffs/{scan}.txt",
    shell:
        "python scripts/process_scores.py {input} > {output}"
