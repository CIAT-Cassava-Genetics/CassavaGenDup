from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import distance
from sklearn import metrics
from sklearn.cluster import DBSCAN

# Fixed analysis parameters
SMISS = 10
EPS = 0.15
MIN_SAMPLES = 2

SNP_FILE = Path("./data/SNPY.csv")
GENOTYPE_FILE = Path("./data/genotypedb.csv")
GENOTYPE_CODES = {"XX": 0, "XY": 1, "YY": 2, "NoCall": 3}


def duplicated_analysis(gt_input_matrix=None, namesdb_input=None):
    """
    
    Identifies duplicate and unique samples in a genotype dataset using
    pairwise Hamming distances and DBSCAN clustering.
        
    Args:
        gt_input_matrix: Genotype matrix from loading_chips. 
        namesdb_input: Sample-to-chip mapping from loading_chips with columns
            'Sample_chip', 'Name' and 'Chip_ID'.

    Returns:
        Tuple of:
            - matrix: Per-sample DataFrame with DBSCAN cluster labels.
              Noise samples (label == -1) are assigned 'unique'.
            - sq_dm: Square-form Hamming distance matrix (np.ndarray).
            - namesdb_tmp: Sample name mapping for all samples included
              in the analysis.
            - freq_samples: Per-sample percentage of each genotype call
              (XX, XY, YY, NoCall) across the selected SNP set.
            - filter_matrix: Numerically encoded genotype matrix restricted
              to samples that passed the missing data filter.

    Raises:
        FileNotFoundError: If SNPY.csv or genotypedb.csv are not found.
        ValueError: If no SNPs with SNPY_PANEL == 'Y' are found, if no samples
            pass the missing data filter or if DBSCAN receives invalid input.
            
    """

    # Load SNP panel
    if not SNP_FILE.is_file():
        raise FileNotFoundError(f"SNP panel file not found: {SNP_FILE.resolve()}")

    snp_df = pd.read_csv(SNP_FILE)
    snp_filter_list = snp_df.loc[snp_df["SNPY_PANEL"] == "Y", "SNPY_ID"].astype(str).tolist()

    if not snp_filter_list:
        raise ValueError(f"No SNPs with SNPY_PANEL == 'Y' found in {SNP_FILE}.")

    # Load historical genotype database
    if not GENOTYPE_FILE.is_file():
        raise FileNotFoundError(f"Genotype database not found: {GENOTYPE_FILE.resolve()}")

    genotypedb = pd.read_csv(GENOTYPE_FILE)
    namesdb = genotypedb[["Sample_chip", "Name", "Chip_ID"]].copy()

    print("Data loaded successfully")

    # Combine incoming chips with historical database (if chips provided)
    if gt_input_matrix is not None and namesdb_input is not None:
        namesdb_tmp = pd.concat([namesdb_input, namesdb], ignore_index=True)
        genotypedb_tmp = pd.concat([gt_input_matrix, genotypedb], ignore_index=True)
    else:
        namesdb_tmp = namesdb.copy()
        genotypedb_tmp = genotypedb.copy()

    # Filter columns to selected SNPs
    available_snps = [s for s in snp_filter_list if s in genotypedb_tmp.columns]
    missing_snps = set(snp_filter_list) - set(available_snps)
    if missing_snps:
        print(f"Warning: {len(missing_snps)} selected SNPs not found in the genotype matrix and will be skipped: {missing_snps}")

    if not available_snps:
        raise ValueError("None of the selected SNPs are present in the genotype matrix.")

    genotypedb_tmp = genotypedb_tmp.set_index("Sample_chip")
    filter_matrix_snps = genotypedb_tmp.loc[:, available_snps]

    # Compute per-sample call frequencies and filter by missing data threshold
    freq_samples = (
        filter_matrix_snps
        .apply(lambda x: x.value_counts(), axis=1)
        .fillna(0)
        .mul(100 / len(available_snps))
    )

    samples_filter_list = freq_samples[freq_samples["NoCall"] <= SMISS].index.tolist()

    if not samples_filter_list:
        raise ValueError(f"No samples passed the missing data filter (NoCall <= {SMISS}%). ")

    print(f"Number of SNPs included: {len(available_snps)}")
    print(f"Number of samples included: {len(samples_filter_list)}")

    # Encode genotypes numerically
    filter_matrix = (
        filter_matrix_snps.loc[samples_filter_list, :]
        .apply(lambda col: col.map(GENOTYPE_CODES))
        .astype(float)
    )

    # Compute pairwise Hamming distance matrix
    print("Building distance matrix...")
    dm = distance.pdist(filter_matrix.values, "hamming")
    sq_dm = distance.squareform(dm)
    print("Distance matrix created")

    # Run DBSCAN clustering
    print("Finding duplicates...")
    db = DBSCAN(eps=EPS, min_samples=MIN_SAMPLES, metric="precomputed").fit(sq_dm)
    labels = db.labels_

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_unique = int((labels == -1).sum())

    # Compute silhouette score
    try:
        metric = round(metrics.silhouette_score(sq_dm, labels, metric="precomputed"), 4)
        print(f"Silhouette score: {metric}")
    except ValueError:
        metric = None
        print("Silhouette score: not applicable. All samples are unique, no clusters detected.")

    # Assign cluster labels
    filter_matrix = filter_matrix.copy()
    filter_matrix["Labels"] = labels
    filter_matrix["sample"] = filter_matrix.index.tolist()

    labels_df = (
        filter_matrix["Labels"]
        .value_counts()
        .reset_index()
    )
    labels_df.columns = ["Labels", "count"]
    labels_df["Cluster"] = range(len(labels_df))

    filter_matrix = filter_matrix.reset_index(drop=True)
    matrix = filter_matrix.merge(labels_df, on="Labels", how="left")
    matrix = matrix.set_index("sample", drop=False)
    matrix = matrix.sort_values("Cluster", ascending=True)
    matrix["Cluster"] = matrix["Cluster"].replace({0: "unique"}).astype(str)
    matrix = matrix.reset_index(drop=True)
    freq_samples = freq_samples.loc[samples_filter_list].reset_index()

    # Report results
    print(f"Number of clusters detected: {n_clusters}")
    print(f"Number of unique samples: {n_unique}")
    print("duplicate detection completed.")

    return matrix, sq_dm, namesdb_tmp, freq_samples, filter_matrix