from importlib.resources import files, as_file
from typing import List, Optional
import pandas as pd

def load_files(
        na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
        ) -> pd.DataFrame:
    
    resource_folder = [f for f in files(__package__).iterdir()]
    out_df = pd.DataFrame()
    for resource in resource_folder:
        with as_file(resource) as file_path:
            if file_path.name.endswith('.csv'):
                tmp = pd.read_csv(file_path, header=0, index_col=False, na_values=na_values)
                out_df = pd.concat([out_df, tmp], ignore_index = True)
    return out_df
