from pathlib import Path
import pandas as pd

DATA_DIR = Path("./data")
RESULTS_DIR = Path("./results")


def curate_db(delete, input, output):
    """
    
    Removes specific samples from a genotype database using a user-defined
    list and saves the result as a new file in ./results/.

    Args:
        delete: CSV filename in ./data/ with one Sample_chip per row
            and no header.
        input: Filename of the existing genotype database in ./data/.
        output: Filename for the curated database saved in ./results/.

    Raises:
        FileNotFoundError: If the input database or delete list are
            not found in ./data/.
        ValueError: If the delete list is empty or if none of the
            samples to delete are found in the database.
            
    """

    input_path = DATA_DIR / input
    delete_path = DATA_DIR / delete

    if not input_path.is_file():
        raise FileNotFoundError(f"Input database not found: {input_path.resolve()}")

    if not delete_path.is_file():
        raise FileNotFoundError(f"Delete list not found: {delete_path.resolve()}")

    # Load existing database
    genotypedb = pd.read_csv(input_path)
    print(f"Database loaded: {len(genotypedb)} samples")

    # Load delete list
    delete_list = pd.read_csv(delete_path, header=None, encoding='utf-8-sig')[0].astype(str).tolist()
    if not delete_list:
        raise ValueError("The delete list is empty.")

    # Validate that samples to delete exist in the database
    existing_samples = set(genotypedb["Sample_chip"].astype(str))
    not_found = set(delete_list) - existing_samples
    if not_found:
        not_found_list = "\n".join(sorted(not_found))
        print(f"Warning: {len(not_found)} samples in the delete list were not found "
              f"in the database and will be ignored:\n{not_found_list}")

    found = [s for s in delete_list if s in existing_samples]
    if not found:
        raise ValueError("None of the samples in the delete list were found in the database.")

    # Remove samples
    curated = genotypedb[~genotypedb["Sample_chip"].astype(str).isin(found)].copy()
    print(f"Samples removed: {len(found)}")
    print(f"Remaining samples: {len(curated)}")

    # Save to results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    outpath = RESULTS_DIR / output
    curated.to_csv(outpath, index=False)
    print(f"Curated database saved at: {outpath}")

    return curated