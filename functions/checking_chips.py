from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap

GENOTYPE_CODES = {"XX": 0, "XY": 1, "YY": 2, "NoCall": 3}
GENOTYPE_LABELS = {0: "XX", 1: "XY", 2: "YY", 3: "NoCall"}
SNP_FILE = Path("./data/SNPY.csv")
RESULTS_DIR = Path("./results")


def checking_chips(gt_input_matrix, chip_ids, figsize=(20, 25)):
    """
    
    Quality control of chip data: computes per-sample genotype-call frequencies 
    and saves a heatmap and frequency table for each chip evaluated.

    Args:
        gt_input_matrix: Genotype matrix from loading_chips. 
        chip_ids: List of chip IDs to evaluate.
        figsize: Figure size in inches for each heatmap. Defaults to (20, 25).

    Returns:
        DataFrame with per-sample genotype call frequencies (%) across all
        evaluated chips combined.

    Raises:
        FileNotFoundError: If ./data/SNPY.csv does not exist.
        ValueError: If a chip_id is not found in gt_input_matrix, or if no
            selected SNPs (SNPY_PANEL == 'Y') are found in the matrix.
    
    """
    # Load and validate SNPY panel file
    if not SNP_FILE.is_file():
        raise FileNotFoundError(f"SNPY panel file not found: {SNP_FILE.resolve()}")

    snp_df = pd.read_csv(SNP_FILE)
    snp_filter_list = snp_df.loc[snp_df["SNPY_PANEL"] == "Y", "SNPY_ID"].astype(str).tolist()

    if not snp_filter_list:
        raise ValueError(f"No SNPs Selected with SNPY_PANEL == 'Y' found in {SNP_FILE}.")

    # Validate that selected SNPs exist in the matrix
    available_snps = [s for s in snp_filter_list if s in gt_input_matrix.columns]
    missing_snps = set(snp_filter_list) - set(available_snps)
    if missing_snps:
        print(f"Warning: {len(missing_snps)} selected SNPs not found in gt_input_matrix and will be skipped: {missing_snps}")

    if not available_snps:
        raise ValueError("None of the selected SNPs are present in gt_input_matrix.")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Color palette
    palette = [sns.color_palette("Set2")[i] for i in [0, 1, 2, 7]]
    cmap = ListedColormap(palette)
    norm = BoundaryNorm([0, 1, 2, 3, 4], ncolors=4)

    freq_list: list[pd.DataFrame] = []

    for chip_id in chip_ids:
        print(f"Checking chip: {chip_id}")

        # Filter matrix by chip
        filtered = gt_input_matrix[gt_input_matrix["Chip_ID"] == chip_id].copy()
        if filtered.empty:
            raise ValueError(f"Chip ID '{chip_id}' not found in gt_input_matrix.")

        # Compute per-sample genotype frequency
        freq = (
            filtered.set_index(["Sample_chip", "Chip_ID"])
            .loc[:, available_snps]
            .apply(lambda x: x.value_counts(), axis=1)
            .fillna(0)
            .mul(100 / len(available_snps))
        )
        freq_list.append(freq)

        # Encode genotypes numerically for heatmap
        encoded = (
            filtered.set_index("Sample_chip")
            .loc[:, available_snps]
            .apply(lambda col: col.map(GENOTYPE_CODES))
            .astype(float)
        )

        # Plot heatmap
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            encoded,
            ax=ax,
            cmap=cmap,
            norm=norm,
            cbar=False,
            xticklabels=encoded.columns.tolist(),
            yticklabels=encoded.index.tolist(),
        )

        patches = [mpatches.Patch(color=palette[i], label=GENOTYPE_LABELS[i]) for i in range(4)]
        ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)

        # Save heatmap
        heatmap_path = RESULTS_DIR / f"{chip_id}_heatmap.png"
        fig.savefig(heatmap_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Heatmap saved at: {heatmap_path}")

    # Concatenate frequency tables from all chips
    freq_all = pd.concat(freq_list)

    # Save frequency table
    if len(chip_ids) == 1:
        csv_filename = f"{chip_ids[0]}_frequency.csv"
    else:
        csv_filename = "multiple_chips_frequency.csv"

    csv_path = RESULTS_DIR / csv_filename
    freq_all.to_csv(csv_path, index=True)
    print(f"Frequency table saved at: {csv_path}")

    return freq_all