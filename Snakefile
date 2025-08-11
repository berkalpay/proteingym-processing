checkpoint align:
    input:
        "data/metadata/{scan}.csv"
    output:
        temp("data/EVcouplings/{scan}/{bitscore}.a2m")
    shell:
        "python scripts/align.py {input} {params.bitscore} > {output}"
    conda:
        "/n/groups/marks/projects/marks_lab_and_oatml/ProteinGym2/EVCouplings/envs/pg2_evc"


def alignment_decision(summary_filepath: str) -> float:
    import pandas as pd

    bitscore = 0  # TODO: implement

    # Extract alignment statistics
    alignment_stats = pd.read_csv(f"{job_name_prefix}_job_statistics_summary.csv")
    assert len(alignment_stats) == 1
    perc_cov = alignment_stats.at[0, "perc_cov"]
    num_seqs = alignment_stats.at[0, "num_seqs"]
    neff_over_l = alignment_stats.at[0, "N_eff"] / alignment_stats.at[0, "seqlen"]

    # Adjust alignment parameters and rerun if necessary
    if (perc_cov < 0.7 or num_seqs < 100 or neff_over_l < 1) and bitscore - 0.05 > 0:
        return bitscore - 0.05
        # TODO: stop if Neff/L < 1 after three tries
    if neff_over_l > 100 and bitscore + 0.05 < 2:
        return bitscore + 0.05  # TODO: skip 1.0
    else:
        return False


rule check_alignment:
    input:
        "data/metadata/{scan}.csv"
    output:
        "data/final_alignments/{scan}.a2m"
    run:
        # TODO: get last alignment
        new_bitscore = alignment_decision()
        if not new_bitscore:
            import shutil
            shutil.copy(aln, output[0])
        else:
            checkpoints.align.get(scan=wildcards.scan, bitscore=new_bitscore)


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
