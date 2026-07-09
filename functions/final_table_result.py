from pathlib import Path
import pandas as pd

PASSPORT_FILE = Path("./data/passportdb.csv")
RESULTS_DIR = Path("./results")


def final_table_result(namesdb_tmp, matrix, freq_samples, filename):
    """
    
    Build the final curated table by merging clustering results, genotype
    call frequencies, sample names, and passport data.
        
    Args:
        namesdb_tmp: Sample name mapping with columns
            'Sample_chip', 'Name' and 'Chip_ID' for all samples.
        matrix: Per-sample DataFrame with columns 'Cluster' and 'count'.
            Contains only samples that passed the missing data filter.
        freq_samples: Per-sample genotype call frequencies. Contains
            only samples that passed the missing data filter.
        filename: Output filename for the Excel file.

    Returns:
        DataFrame with the final table.

    Raises:
        FileNotFoundError: If passportdb.csv is not found.
        
    """

    # Load passport database
    if not PASSPORT_FILE.is_file():
        raise FileNotFoundError(f"Passport database not found: {PASSPORT_FILE.resolve()}")

    passportdb = pd.read_csv(PASSPORT_FILE)

    # Build base table
    
    matrix["Cluster"] = matrix["Cluster"].astype(str)
    
    clusterdb = (
        namesdb_tmp
        .merge(matrix[["sample", "Cluster", "count"]], left_on="Sample_chip", right_on="sample", how="inner")
        .drop(columns=["sample"])
        .merge(freq_samples, on="Sample_chip", how="left")
    )

    # Merge passport data
    tmp_data = (
        clusterdb
        .merge(passportdb, left_on="Sample_chip", right_on="Name_Chip", how="left", indicator=True)
        .rename(columns={"_merge": "passport_status"})
    )
    tmp_data["passport_status"] = tmp_data["passport_status"].map(
        {"left_only": "No_passport_data", "both": "Ready"}
    )

    # Drop redundant passport columns
    cols_to_drop = [c for c in ["Name_Chip", "passport_name"] if c in tmp_data.columns]
    final_table = tmp_data.drop(columns=cols_to_drop)

    # Save final table as Excel
    
    id_cols = ["Sample_chip", "Name", "Chip_ID", "Cluster", "count"]
    freq_cols = ["NoCall", "XX", "XY", "YY"]
    passport_cols = [c for c in final_table.columns if c not in id_cols + freq_cols]
    final_table = final_table[id_cols + freq_cols + passport_cols]
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    outpath = RESULTS_DIR / filename
    final_table.to_excel(outpath, index=False)
    print(f"Final table saved at: {outpath}")

    return final_table