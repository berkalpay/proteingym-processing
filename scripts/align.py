import os
import sys
import pandas as pd


# Extract relevant metadata
metadata_path = sys.argv[1]
metadata = pd.read_csv(metadata_path)
dms_id = metadata.at[0, "dms_id"]
virus = metadata.at[0, "taxon"] == "Virus"
target_seq = metadata.at[0, "target_seq"]

# Set parameters and run alignment
theta = 1 - (0.01 if virus else 0.2)
bitscore = 0.2 if virus else 0.7
target_seq_path = f"temp/{dms_id}"
with open(target_seq_path, "w") as f:
    f.write(">{dms_id}\n")
    f.write(target_seq)
job_name_prefix = f"data/{dms_id}/bit_{bitscore}_theta_{theta}_colcov_70"
os.system(
    f"evcouplings --protein {dms_id} -b {bitscore} -s {target_seq_path}"
    f"--theta {theta} --colcov 70 --prefix {job_name_prefix} configs/align.txt --yolo"
)
# TODO: use --stages arg instead of separate configs?

# Extract alignment statistics
alignment_stats = pd.read_csv(f"{job_name_prefix}_job_statistics_summary.csv")
assert len(alignment_stats) == 1
perc_cov = alignment_stats.at[0, "perc_cov"]
num_seqs = alignment_stats.at[0, "num_seqs"]
neff_over_l = alignment_stats.at[0, "N_eff"] / alignment_stats.at[0, "seqlen"]

# Adjust alignment parameters and rerun if necessary
if (perc_cov < 0.7 or num_seqs < 100 or neff_over_l < 1) and bitscore - 0.05 > 0:
    pass  # TODO: decrease bitscore by 0.05 and rerun
    # TODO: if Neff/L < 1 after three tries, take the alignment with the highest Neff/L
if neff_over_l > 100 and bitscore + 0.05 < 2:
    pass  # TODO: increase the bitscore by 0.05 (skipping bitscore=1.0) and rerun
