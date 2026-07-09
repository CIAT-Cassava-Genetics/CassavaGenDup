from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap

SNP_FILE = Path("./data/SNPY.csv")
RESULTS_DIR = Path("./results")
GENOTYPE_LABELS = {0: "XX", 1: "XY", 2: "YY", 3: "NoCall"}


def heatmap_samples(sample_list_file, matrix, filename, dpi, long):
    """
    
    Generates a genotype heatmap for a user-defined list of samples.
        
    Args:
        sample_list_file: CSV filename located in ./data/ with one
            Sample_chip per row and no header.
        matrix: Per-sample DataFrame from duplicated_analysis with column
            'sample' and SNP genotype codes.
        filename: Output filename.
        dpi: Image resolution in dots per inch.
        long: Figure height in inches.

    Raises:
        FileNotFoundError: If SNPY.csv or sample_list_file are not found.
        ValueError: If sample_list is empty or no selected SNPs are present
            in the matrix.
            
    """

    # Validate input files
    if not SNP_FILE.is_file():
        raise FileNotFoundError(f"SNP panel file not found: {SNP_FILE.resolve()}")

    # Load selected SNPs
    snp_df = pd.read_csv(SNP_FILE)
    snp_filter_list = snp_df.loc[snp_df["SNPY_PANEL"] == "Y", "SNPY_ID"].astype(str).tolist()
    snp_cols = [c for c in snp_filter_list if c in matrix.columns]

    if not snp_cols:
        raise ValueError("None of the selected SNPs are present in the matrix.")

    # Load sample list 
    sample_file = Path("./data") / sample_list_file
    if not sample_file.is_file():
        raise FileNotFoundError(f"Sample list file not found: {sample_file.resolve()}")

    sample_list = pd.read_csv(sample_file, header=None, encoding='utf-8-sig')[0].astype(str).tolist()
    if not sample_list:
        raise ValueError("The sample list is empty.")

    # Filter matrix to requested samples
    sample_snpy = (
        matrix.set_index("sample")
        .reindex(sample_list)
        [snp_cols]
        .apply(pd.to_numeric, errors="coerce")
    )

    # Color palette
    palette = [sns.color_palette("Set2")[i] for i in [0, 1, 2, 7]]
    cmap = ListedColormap(palette)
    norm = BoundaryNorm(boundaries=[0, 1, 2, 3, 4], ncolors=4)

    # Plot heatmap
    fig, ax = plt.subplots(figsize=(18, long))
    sns.heatmap(
        sample_snpy,
        ax=ax,
        cmap=cmap,
        norm=norm,
        cbar=False,
        xticklabels=sample_snpy.columns.tolist(),
        yticklabels=sample_snpy.index.tolist(),
    )

    # Legend
    patches = [mpatches.Patch(color=palette[i], label=GENOTYPE_LABELS[i]) for i in range(4)]
    ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)

    # Save heatmap
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    outpath = RESULTS_DIR / filename
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Samples heatmap saved at: {outpath}")