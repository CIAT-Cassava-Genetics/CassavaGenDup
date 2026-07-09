# CassavaGenDup

**CassavaGenDup** is a Python/Jupyter Notebook tool developed for cassava varietal identification and duplicate detection using SNP genotyping data from the SNPY-Chip.

The pipeline allows to load genotype files, evaluate genotype quality, identify duplicated and unique samples, generate final reports, visualize genotype profiles, export VCF files, and perform kinship-based relationship analysis.

This repository includes an initial curated reference database of approximately 12,000 cassava samples.

## How to use

### 1. Clone or download this repository

```bash
git clone https://github.com/yourusername/CassavaGenDup.git
cd CassavaGenDup
```

### 2. Create a virtual environment

Creating a virtual environment is recommended to avoid conflicts with other Python installations.

CassavaGenDup is recommended to run with Python 3.11. 

```bash
conda create -n duplicated_analysis_env python=3.11
conda activate duplicated_analysis_env
```

### 3. Install the required Python packages

```bash
pip install -r requirements.txt
```

### 4. Add input files

Place the genotype input files to be analyzed in the `input/` folder.

### 5. Open the notebook

Open JupyterLab from your local computer or  server, navigate to the `CassavaGenDup` folder, and open:

```text
CassavaGenDup.ipynb
```

### 6. Run the analysis

Run the notebook sections step by step.

## Input files

Input genotype files should be placed in the `input/` folder.

## Reference database

The `data/` folder contains the curated reference files used by the pipeline:

* `genotypedb.csv`: curated genotype database.
* `passportdb.csv`: curated passport database.
* `SNPY.csv`: SNP marker information.

The initial release includes approximately 12,000 cassava samples.

## Main outputs

The pipeline can generate:

* Genotype quality summaries.
* Duplicate and unique sample classification.
* Final Excel reports.
* Genotype heatmaps.
* VCF files.
* Kinship heatmaps and pairwise relationship tables.

Output files are saved in the `results/` folder.

## Citation

If you use CassavaGenDup or the SNPY-Chip data structure in your research, please cite the associated manuscript:

Gomez-Martinez L., Barrera V., Bohórquez-Chaux A., Ovalle T., Perea C., Leiva A.M., Narvaez M., Vélez M., Labarta R., Peña-Venegas C.P., Ocampo J., Dufour D., Newby J., Gimode W., and Becerra López-Lavalle L.A.
**A global genomic map of cassava diversity: development of the SNPY-Chip and the CassavaGenDup tool for varietal identification.**
Manuscript in preparation.

## Authors

CassavaGenDup was developed as part of cassava genetic diversity and varietal identification research at Alliance Bioversity International and CIAT.

For questions about this tool, please contact the Cassava Genetics Laboratory team at Alliance Bioversity International and CIAT: [alliance-cassavagenetics@cgiar.org]

## License

This project is licensed under the GNU General Public License v3.0 (GPLv3). See the [LICENSE](LICENSE) file for details.