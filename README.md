# code_scbasset

# Source Folder

The source folder contains a variety of Jupyter Notebooks as well as Python and Bash scripts designed for the preprocessing and handling of genomic data. Below is a detailed description of each file present in this directory:

## Preprocessing

### CAGE_preprocessing
- **Description**: Converts an original CSV file (containing enhancer and promoter data) into a BED file listing the peak regions (promoter and enhancer). Also generates an H5 file containing the count matrix.
- **Input**: CSV file (enhancer and promoter).
- **Output**:
  - BED file of peaks.
  - H5 file (count matrix).

### Concat_enhancer_promoteur
- **Description**: Concatenates two CSV files, typically those of promoters and enhancers.
- **Input**: Two CSV files (promoter and enhancer).
- **Output**: A combined CSV file.

### Controle_shuffle
- **Description**: Randomly shuffles the values in the original CSV file to create a control.
- **Input**: Original CSV file.
- **Output**: Shuffled CSV file.

### SCAFE_preprocessing
- **Description**: From a sparse matrix and two text files containing column and row information, creates an H5 file containing the count matrix.
- **Input**:
  - Sparse matrix.
  - Text files for columns and rows.
- **Output**: H5 file (count matrix).

### Testing-random
- **Description**: A no-op script intended for various tests.
- **Function**: None (testing purposes only).

### make_anndata
- **Description**: Generates the final file necessary for creating training files. Allows filtering of the number of enhancers, promoters, etc.
- **Function**: Data preparation for training.
