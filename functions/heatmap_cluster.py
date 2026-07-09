from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import BoundaryNorm, ListedColormap
 
SNP_FILE = Path("./data/SNPY.csv")
RESULTS_DIR = Path("./results")
GENOTYPE_LABELS = {0: "XX", 1: "XY", 2: "YY", 3: "NoCall"}
 
 
def heatmap_cluster(cluster_list, matrix, filename, dpi, figsize=(20, 10)):
    """
    
    Generates a genotype heatmap for one or more specified clusters from the
    duplicate analysis matrix.
    
    Args:
        cluster_list: List of cluster IDs to include.
        matrix: Per-sample DataFrame from duplicated analysis with columns
            'Cluster' and SNP genotype codes.
        filename: Output filename.
        dpi: Image resolution in dots per inch.
        figsize: Figure size in inches. Defaults to (20, 10).
    
    Raises:
        FileNotFoundError: If SNPY.csv is not found.
        ValueError: If cluster_list is empty or contains IDs not present
            in the matrix.
            
    """
    
    # Validate input files
    if not SNP_FILE.is_file():
       raise FileNotFoundError(f"SNP panel file not found: {SNP_FILE.resolve()}")    
    
    # Validate requested clusters exist in matrix
    available_clusters = set(matrix["Cluster"].unique())
    missing_clusters = [c for c in cluster_list if str(c) not in available_clusters]
    if missing_clusters:
        raise ValueError(f"The following clusters are not present in the matrix: {missing_clusters}")
    
    # Load selected SNPs
    snp_df = pd.read_csv(SNP_FILE)
    snp_filter_list = snp_df.loc[snp_df["SNPY_PANEL"] == "Y", "SNPY_ID"].astype(str).tolist()
    snp_cols = [c for c in snp_filter_list if c in matrix.columns]
    
    if not snp_cols:
        raise ValueError("None of the selected SNPs are present in the matrix.")
 
    # Filter matrix to requested clusters
    subset = matrix[matrix["Cluster"].isin([str(c) for c in cluster_list])].copy()
    cluster_values = subset.set_index("sample")[snp_cols].apply(pd.to_numeric, errors="coerce")
 
    # Color palette
    palette = [sns.color_palette("Set2")[i] for i in [0, 1, 2, 7]]
    cmap = ListedColormap(palette)
    norm = BoundaryNorm(boundaries=[0, 1, 2, 3, 4], ncolors=4)
 
    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        cluster_values,
        ax=ax,
        cmap=cmap,
        norm=norm,
        cbar=False,
        xticklabels=cluster_values.columns.tolist(),
        yticklabels=cluster_values.index.tolist(),
    )
 
    # Legend
    patches = [mpatches.Patch(color=palette[i], label=GENOTYPE_LABELS[i]) for i in range(4)]
    ax.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
 
    # Save heatmap
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    outpath = RESULTS_DIR / filename
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Cluster heatmap saved at: {outpath}")