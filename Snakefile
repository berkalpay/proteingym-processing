checkpoint align:
    input:
        "data/metadata/{scan}.csv",
    output:
        temp("data/EVcouplings/{scan}/{bitscore}.a2m"),
    conda:
        "/n/groups/marks/projects/marks_lab_and_oatml/ProteinGym2/EVCouplings/envs/pg2_evc"
    shell:
        "python scripts/align.py {input} {params.bitscore} > {output}"


def alignment_decision(summary_filepath: str) -> float:
    import pandas as pd

    bitscore = float(
        summary_filepath.split(os.sep)[-1].split("_alignment_statistics.csv")[0]
    )
    STEP = 0.05

    # Extract alignment statistics
    alignment_stats = pd.read_csv(f"{job_name_prefix}_job_statistics_summary.csv")
    assert len(alignment_stats) == 1
    perc_cov = alignment_stats.at[0, "perc_cov"]
    num_seqs = alignment_stats.at[0, "num_seqs"]
    neff_over_l = alignment_stats.at[0, "N_eff"] / alignment_stats.at[0, "seqlen"]

    # Adjust alignment parameters and rerun if necessary
    if (perc_cov < 0.7 or num_seqs < 100 or neff_over_l < 1) and bitscore - STEP > 0:
        return bitscore - STEP
    if neff_over_l > 100 and bitscore + STEP < 2:
        new_bitscore = bitscore + 0.05
        if new_bitscore == 1:
            return bitscore + 2 * STEP
        return new_bitscore
    else:
        return False


rule check_alignment:
    input:
        "data/metadata/{scan}.csv",
    output:
        "data/final_alignments/{scan}.a2m",
    run:
        alignments_dir = "data/EVcouplings/{wildcards.scan}/"
        files_to_times = {}

        for bitscore in os.listdir(alignments_dir):
            bitscore_dir = os.path.join(alignments_dir, bitscore)
            for filename in [
                f
                for f in os.listdir(bitscore_dir)
                if f.endswith("_alignment_statistics.csv")
            ]:
                file = os.path.join(bitscore_dir, filename)
                files_to_times[file] = os.stat(file).st_mtime
        latest_file = max(files_to_times, key=lambda x: files_to_times[x])

        new_bitscore = alignment_decision(latest_file)
        if not new_bitscore:
            import shutil

            shutil.copy(
                latest_file.replace("_alignment_statistics.csv", ".a2m"), output[0]
            )
        else:
            checkpoints.align.get(scan=wildcards.scan, bitscore=new_bitscore)


rule process_scores:
    input:
        "data/metadata/{scan}.csv" "data/raw_scores/{scan}.csv",
    output:
        "data/scores/{scan}.csv" "data/binarization_cutoffs/{scan}.txt",
    shell:
        "python scripts/process_scores.py {input} > {output}"
