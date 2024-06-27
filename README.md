# Quick start





# code_scbasset

# Architecture

```
├── README.md
├── source
│   ├── batch_script
│   │   ├── make_train_data.sh
│   │   ├── modisco_lite.sh
│   │   ├── shap.sh
│   │   └── train_model.sh
│   ├── DeepSHAP
│   │   ├── shap_functions.py
│   │   └── shap_multiprocessing.py
│   ├── scBasset
│   │   ├── analyse
│   │   │   ├── analyse_function.py
│   │   │   ├── AUC_Loss.py
│   │   │   ├── embedding.py
│   │   │   ├── overlapping.py
│   │   └── preprocessing
│   │       ├── CAGE_preprocessing.py
│   │       ├── Controle_shuffle.py
│   │       ├── preprocessing_function.py
└── test
```


# scBasset (source folder)

## Preprocessing Folder

This folder contains a number of Python scripts for pre-processing data in order to train the model. Below you will find descriptions of each script:

### CAGE_preprocessing.py
- **Description**: Converts an original CSV file (row: DNA position, column: sample), containing enhancer and promoter data, into a BED file listing the promoter or enhancer regions. Also generates an H5 file containing the count matrix and another H5 file with the count matrix filtered (rate to be defined).
- **Input**: CSV file (enhancer and promoter).
- **Output**:
  - BED file of peaks.
  - H5 file (count matrix).
  - H5 file (count matrix filtered).

### Controle_shuffle.ipynb
- **Description**: Randomly shuffles the values in the original CSV file to create a control.
- **Input**: Original CSV file.
- **Output**: Shuffled CSV file.

## Analysis Folder

The analysis directory includes several Python scripts designed for visualizing various aspects of model training and evaluation. Below you will find descriptions of each script:

### AUC_Loss.py
- **Description**: Allows visualization of the loss and AUC (Area Under the Curve) for both training and validation phases of each model training session.
- **Input**: history.pickle file
- **Output**: PNG with two plots (AUC and Loss)

### Embedding.py
- **Description**: Script to perform an embedding of the weights of the last dense layer of the model. Dimensionality reduction is then performed, allowing visualization of the different samples colored by (patient, cell lines, or hormonal status) in different plots (PNG).
- **Input**: Count matrix (H5 file), trained model, CSV with hormone status information
- **Output**: 3 embedding PNG plots (leiden, patients vs cell lines, hormone status)

### Overlapping.py
- **Description**: Script that calculates the percentage of overlap between chromosome positions (extended to 1344 bp) within a BED file. 
Possibility to visualize the global overlap distribution through a histogram as well as a chromosome-by-chromosome visualization. 
- **Input**: BED file
- **Output**: Overlapping percentage, two PNG plots


# DeepSHAP (source folder)

### shap_multiprocessing.py
- **Description**: Script that calculates the importance of each nucleotide (Shapley values) for the final model prediction for a given number of input sequences 
- **Input**: Sequence to analyze, trained model, number of cores for parallelization 
- **Output**: Two NPZ files (one with original sequences, one with the corresponding Shapley values)




