import os
import sys
import pandas as pd

metadata_path = sys.argv[1]
extracted_metadata_path = sys.argv[2]

metadata = pd.read_csv(metadata_path)
for dms_id in metadata["DMS_id"]:
    extracted = metadata.loc[metadata["DMS_id"] == dms_id]
    if not (
        os.path.exists(extracted_metadata_path)
        and pd.read_csv(extracted_metadata_path) == extracted
    ):
        extracted.to_csv(extracted_metadata_path + f"{dms_id}.csv", index=False)
