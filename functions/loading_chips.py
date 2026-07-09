from pathlib import Path 
import numpy as np
import pandas as pd

valid_calls = {"XX", "XY", "YY", "NoCall"}

replace_val = {"Invalid": "NoCall", "No Call": "NoCall",}

drop_val = {"NTC", "Empty"}


def loading_chips(final_call, assay_name, sample_name, skiprows):    
    """
    
    Loads and pre-processes SNP genotyping files (CSV format) from the ./input/
    directory and returns the data structures required by the CassavaGenDup pipeline.
    
    Args:
    final_call: Column name for the genotype call.
    assay_name: Column name for the SNP assay identifier.
    sample_name: Column name for the sample name.
    skiprows: Header rows to skip before column names.

    Returns:
    Tuple of (chip_all, gt_input_matrix, namesdb_input).

    Raises:
    FileNotFoundError: If ./input/ does not exist.
    ValueError: If no CSVs found or duplicate (Sample_chip, assay) pairs detected.
    
    """
   
    input_dir = Path("./input")

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir.resolve()}")

    #Load SNP genotyping files
    csv_files = sorted(input_dir.glob("*.csv"))
    if not csv_files:
        raise ValueError(f"No CSV files found in {input_dir.resolve()}")

    print(f"Files detected: {len(csv_files)}")

    frames: list[pd.DataFrame] = []

    for csv_path in csv_files:
        chip_id = csv_path.stem
        print(f"Loading chip: {chip_id}")

        raw = pd.read_csv(
            csv_path,
            skiprows=skiprows,
            header=0,
            dtype=object,
        )

        # Standardise genotype calls and remove uninformative rows
        raw[final_call] = raw[final_call].replace(replace_val)
        raw.loc[raw[sample_name].isin(drop_val), sample_name] = np.nan

        # Build unique sample identifier
        raw["Sample_chip"] = raw[sample_name].astype(str) + "_" + chip_id
        raw["Chip_ID"] = chip_id

        # Keep only the columns needed and drop rows with any missing value
        subset = (
            raw[[assay_name, "Sample_chip", final_call, sample_name, "Chip_ID"]]
            .dropna()
        )
        
        # Validate that genotype calls are expected values
        unexpected = set(subset[final_call].unique()) - valid_calls
        if unexpected:
            raise ValueError(
                f"Chip '{chip_id}' contains unexpected genotype values in "
                f"column '{final_call}': {unexpected}. "
                f"Expected values are: {valid_calls}."
            )

        frames.append(subset)

    chip_all = pd.concat(frames, ignore_index=True)

    # Detect duplicate (Sample_chip, assay) pairs
    duplicate_mask = chip_all.duplicated(subset=["Sample_chip", assay_name], keep=False)
    if duplicate_mask.any():
        duplicates = (
            chip_all.loc[duplicate_mask, ["Sample_chip", assay_name]]
            .drop_duplicates()
            .head(10)
            .to_string(index=False)
        )
        raise ValueError(
            f"Duplicate (Sample_chip, assay) pairs detected:\n{duplicates}"
        )

    # Pivot to sample × SNP matrix
    gt_input_matrix = chip_all.pivot(
        index=["Sample_chip", sample_name, "Chip_ID"],
        columns=assay_name,
        values=final_call,
    ).reset_index()
    
    gt_input_matrix.columns.name = None

    #Extract simple mapping table
    namesdb_input = gt_input_matrix[["Sample_chip", sample_name, "Chip_ID"]].copy()

    print(f"Loaded {len(namesdb_input)} samples across {len(csv_files)} chips.")

    return chip_all, gt_input_matrix, namesdb_input