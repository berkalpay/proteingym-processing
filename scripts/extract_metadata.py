import sys
import pandas as pd

metadata_path = sys.argv[1]
dms_id = sys.argv[2]
extracted_metadata_path = sys.argv[3]

metadata = pd.read_csv(metadata_path)
extracted = metadata.loc[metadata["DMS_id"] == dms_id]
if len(extracted) == 0:
    raise ValueError("DMS ID {dms_id} isn't in the given metadata file")
elif len(extracted) > 1:
    raise ValueError("DMS ID {dms_id} is in the given metadata file multiple times")
extracted.to_csv(extracted_metadata_path, index=False)
