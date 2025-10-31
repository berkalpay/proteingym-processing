import os
import sys
import pandas as pd


metadata_path = sys.argv[1]
bitscore = sys.argv[2]

# Extract relevant metadata
metadata = pd.read_csv(metadata_path)
dms_id = metadata.at[0, "DMS_id"]
virus = metadata.at[0, "taxon"] == "Virus"
target_seq = metadata.at[0, "target_aa_seq"]

# Run alignment
theta = 1 - (0.01 if virus else 0.2)
os.makedirs("data/temp/", exist_ok=True)
target_seq_path = f"data/temp/{dms_id}"
with open(target_seq_path, "w") as f:
    f.write(">{dms_id}\n")
    f.write(target_seq)
job_name_prefix = f"data/EVcouplings/{dms_id}/{bitscore}"
os.system(
    f"evcouplings --protein {dms_id} -b {bitscore} -s {target_seq_path}"
    f"--theta {theta} --colcov 70 --prefix {job_name_prefix} configs/align.txt --yolo"
)  # TODO: use --stages arg instead of separate configs?
