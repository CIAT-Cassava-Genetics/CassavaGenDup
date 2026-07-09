from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

RESULTS_DIR = Path("./results")


def relationship(phi):

    if phi < 0.0442:
        return "Unrelated"
    if phi < 0.0884:
        return "3rd Degree"
    if phi < 0.1768:
        return "2nd Degree"
    if phi < 0.3536:
        return "1st Degree"
    return "Twin"


def kinship_plot(kinship_file, kinship_output, matrix, dpi=600):
    """
    
    Generates a kinship heatmap and pairwise relationship table from VCFtools
    relatedness2 output.
        
    Args:
        kinship_file: Input filename from ./data/
        kinship_output: Output CSV filename for the pairwise table. 
            The heatmap is saved with the same name plus '_heatmap.tiff'.
        matrix: Per-sample DataFrame from duplicated_analysis with
            columns 'sample' and 'Cluster'.
        dpi: Image resolution in dots per inch.

    Raises:
        FileNotFoundError: If kinship_file is not found in ./data/.
        
    """
 
    in_path = Path("./data") / kinship_file
    if not in_path.is_file():
        raise FileNotFoundError(f"Kinship file not found: {in_path.resolve()}")
 
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    heatmap_path = RESULTS_DIR / (Path(kinship_output).stem + "_heatmap.tiff")
    out_csv = RESULTS_DIR / kinship_output
 
    # Load and process kinship data
    kinship_long = pd.read_csv(in_path, sep="\t")
    kinship_long["PHI"] = kinship_long["RELATEDNESS_PHI"].clip(lower=0)
    kinship_long = kinship_long[["INDV1", "INDV2", "PHI"]]
 
    # Pivot to square matrix
    kinship = kinship_long.pivot(index="INDV1", columns="INDV2", values="PHI")
 
    # Generate clustermap
    cmap = sns.cubehelix_palette(start=0.5, rot=-0.5, as_cmap=True)
 
    cluster_fig = sns.clustermap(
        kinship,
        cmap=cmap,
        figsize=(20, 20),
        dendrogram_ratio=(0.1, 0.1),
        cbar_pos=(1.02, 0.3, 0.02, 0.4),
    )
 
    cluster_fig.fig.savefig(heatmap_path, dpi=dpi, bbox_inches="tight")
    plt.close(cluster_fig.fig)
    print(f"Kinship heatmap saved at: {heatmap_path}")
 
    # Build pairwise table (upper triangle only, excluding self-comparisons)
    pairwise = (
        kinship
        .where(np.triu(np.ones(kinship.shape), k=1).astype(bool))
        .stack()
        .reset_index()
    )
    pairwise.columns = ["INDV1", "INDV2", "PHI"]
    pairwise = pairwise.sort_values("PHI", ascending=True)
 
    # Merge cluster labels for both samples
    cluster_map = matrix.set_index("sample")["Cluster"]
    pairwise["INDV1_label"] = pairwise["INDV1"].map(cluster_map)
    pairwise["INDV2_label"] = pairwise["INDV2"].map(cluster_map)
 
    # Classify relationship
    pairwise["relationship"] = pairwise["PHI"].apply(relationship)
 
    pairwise.to_csv(out_csv, index=False)
    print(f"Kinship table saved at: {out_csv}")
 
    return pairwise