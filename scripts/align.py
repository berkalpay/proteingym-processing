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
job_name_prefix = f"data/EVcouplings/{dms_id}/bit_{bitscore}_theta_{theta}_colcov_70"
os.system(
    f"evcouplings --protein {dms_id} -b {bitscore} -s {target_seq_path}"
    f"--theta {theta} --colcov 70 --prefix {job_name_prefix} configs/align.txt --yolo"
)
# TODO: use --stages arg instead of separate configs?
