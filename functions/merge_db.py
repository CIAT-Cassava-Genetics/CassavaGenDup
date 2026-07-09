from pathlib import Path
import pandas as pd

DATA_DIR = Path("./data")
RESULTS_DIR = Path("./results")


def merge_db(gt_input_matrix, namesdb_input, input, output):
    """

    Merges new chip data with an existing genotype database and saves the
    result as a new file in ./results/.
        
    Args:
        gt_input_matrix: Genotype matrix from loading_chips.
        namesdb_input: Sample name mapping from loading_chips with columns
            'Sample_chip', 'Name' and 'Chip_ID'.
        input: Filename of the existing genotype database in ./data/.
        output: Filename for the merged database saved in ./results/.

    Raises:
        FileNotFoundError: If the input database file is not found in ./data/.
        ValueError: If duplicate Sample_chip are found between the new
            chips and the existing database.
            
    """

    input_path = DATA_DIR / input
    if not input_path.is_file():
        raise FileNotFoundError(f"Input database not found: {input_path.resolve()}")

    # Load existing database
    genotypedb = pd.read_csv(input_path)
    print(f"Existing database loaded: {len(genotypedb)} samples")

    # Check for duplicates
    existing_samples = set(genotypedb["Sample_chip"].astype(str))
    new_samples = set(gt_input_matrix["Sample_chip"].astype(str))
    duplicates = existing_samples & new_samples

    if duplicates:
        duplicate_list = "\n".join(sorted(duplicates))
        raise ValueError(
            f"The following {len(duplicates)} Sample_chip already exist in "
            f"the database and cannot be merged:\n{duplicate_list}"
        )

    # Merge new chips with existing database
    merged = pd.concat([gt_input_matrix, genotypedb], ignore_index=True)
    print(f"New samples added: {len(gt_input_matrix)}")
    print(f"Total samples in merged database: {len(merged)}")

    # Save to results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    outpath = RESULTS_DIR / output
    merged.to_csv(outpath, index=False)
    print(f"Merged database saved at: {outpath}")

    return merged