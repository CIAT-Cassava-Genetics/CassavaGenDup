# Input files

Place the genotype input files to be analyzed by CassavaGenDup in this folder.

The pipeline expects genotype files in `.csv` format. For standard Fluidigm chip files, the notebook uses the following default parameters:

- `final_call = "Final"`
- `assay_name = "Assay"`
- `sample_name = "Name"`
- `skiprows = 15`

Each input file should correspond to one genotyping chip or plate.

