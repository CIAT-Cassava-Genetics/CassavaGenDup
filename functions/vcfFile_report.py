from pathlib import Path
import pandas as pd

SNP_FILE = Path("./data/SNPY.csv")
RESULTS_DIR = Path("./results")
GT_MAP = {0: "0/0", 1: "0/1", 2: "1/1", 3: "./."}


def vcfFile_report(sample_list_file, matrix, output_vcf):
    """
    
    Generates a VCF file for a user-defined list of samples using the
    genotype matrix and SNP metadata from SNPY.csv.
        
    Args:
        sample_list_file: CSV filename in ./data/ with one Sample_chip
            per row and no header.
        matrix: Per-sample DataFrame from duplicated_analysis with column
            'sample' and SNP genotype codes.
        output_vcf: Output VCF filename.

    Raises:
        FileNotFoundError: If SNPY.csv or sample_list_file are not found.
        ValueError: If sample list is empty or no selected SNPs are present
            in the matrix.
     
    """

    # Validate input files
    if not SNP_FILE.is_file():
        raise FileNotFoundError(f"SNP panel file not found: {SNP_FILE.resolve()}")

    sample_file = Path("./data") / sample_list_file
    if not sample_file.is_file():
        raise FileNotFoundError(f"Sample list file not found: {sample_file.resolve()}")

    # Load SNP metadata — filtered to selected SNPs (SNPY_PANEL == 'Y')
    snp_df = pd.read_csv(SNP_FILE)
    snp_filter_list = snp_df.loc[snp_df["SNPY_PANEL"] == "Y", "SNPY_ID"].astype(str).tolist()

    # Load sample list
    sample_list = pd.read_csv(sample_file, header=None, encoding='utf-8-sig')[0].astype(str).tolist()
    if not sample_list:
        raise ValueError("The sample list is empty.")

    # Filter matrix to requested samples and available SNPs
    snp_cols = [c for c in snp_filter_list if c in matrix.columns]
    if not snp_cols:
        raise ValueError("None of the selected SNPs are present in the matrix.")

    report = (
        matrix.set_index("sample")
        .reindex(sample_list)
        [snp_cols]
    )

    # Convert numeric genotype codes to VCF format
    gt_data = report.apply(lambda col: col.map(GT_MAP))

    # Transpose for VCF format 
    gt_data = gt_data.T

    # Build VCF metadata columns from SNPY.csv
    meta = snp_df[snp_df["SNPY_ID"].isin(snp_cols)].copy()
    meta = meta.rename(columns={"SNPY_ID": "ID", "CHROM": "#CHROM", "X": "REF", "Y": "ALT"})
    meta["QUAL"] = 255
    meta["FILTER"] = "PASS"
    meta["INFO"] = "."
    meta["FORMAT"] = "GT"

    # Arrange VCF columns in standard order
    vcf_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    meta = meta[vcf_cols]

    # Merge metadata with genotype data
    vcf_body = meta.merge(
        gt_data,
        left_on="ID",
        right_index=True,
        how="inner",
    )

    # Write VCF file

    header_meta = (
        "##fileformat=VCFv4.2\n"
        "##source=CassavaGenDupPipeline\n"
        "##reference=Manihot_esculenta_v8\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    )
    header_cols = "\t".join(vcf_cols + sample_list) + "\n"

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    outpath = RESULTS_DIR / output_vcf
    
    with open(outpath, "w") as vcf:
        vcf.write(header_meta)
        vcf.write(header_cols)

    vcf_body.to_csv(outpath, sep="\t", mode="a", index=False, header=False)
    print(f"VCF file saved at: {outpath}")