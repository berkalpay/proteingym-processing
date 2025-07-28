import os
import sys
import pandas as pd


metadata_path = sys.argv[1]
metadata = pd.read_csv(metadata_path)

dms_id = metadata.at[0, "dms_id"]
virus = metadata.at[0, "taxon"] == "Virus"
theta = 1 - (0.01 if virus else 0.2)
bitscore = 0.2 if virus else 0.7
target_seq_path = f"temp/{dms_id}"
with open(target_seq_path, "w") as f:
    f.write(">{dms_id}\n")
    f.write(metadata.at[0, "target_seq"])

job_name_prefix = f"proteingym_{dms_id}_bit_{bitscore}_theta_{theta}_colcov_70"
os.system(
    f"evcouplings --protein {dms_id} -b {bitscore} -s {target_seq_path}"
    f"--theta {theta} --colcov 70 --prefix {job_name_prefix} configs/align.txt --yolo"
)
# TODO: use --stages arg instead of separate configs?

# # OUTPUT
# 1. .a2m file
# 2. aln_stats file
# 3. weights file
# 4. align_config file

# # check whether the alignment is good enough
# Check the f: {job_name_prefix}_job_statistics_summary.csv file
# 1. If <70% of the columns covered by at least 70% of the sequences in the final alignment (perc_cov), or num_seqs < 100:
# lower the bitscore by 0.5 (not lower than 0)
# 2. If Neff/L > 100 (Neff divided by seqlen):
# increase the bitscore by 0.5 (not higher than 2, skip 1.0)
# 3. If Neff/L < 1:
# lower the bitscore by 0.5 (not lower than 0)
# if still < 1 after three tries, take the alignment with the highest Neff/L

# For design sequences, as long as the Neff/L is close to 100, it's ok.

# # iterate above till good
