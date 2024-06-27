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

# Source Folder



## Preprocessing Folder

This folder contains a number of python scripts for pre-processing data in order to train the model. Below you will find descriptions of each script:

### CAGE_preprocessing.py
- **Description**: Converts an original CSV file (row: DNA position, column: sample), containing enhancer and promoter data, into a BED file listing the promoter or enhancer regions. Also generates an H5 file containing the count matrix and an other H5 file with the coutn matrx filtered (rate to defined).
- **Input**: CSV file (enhancer and promoter).
- **Output**:
  - BED file of peaks.
  - H5 file (count matrix).
  - H5 file (coutn matrix filtered).

### Controle_shuffle.ipynb
- **Description**: Randomly shuffles the values in the original CSV file to create a control.
- **Input**: Original CSV file.
- **Output**: Shuffled CSV file.

## Analysis Folder

The analysis directory includes several pythons scripts designed for visualizing various aspects of model training and evaluation. Below you will find descriptions of each script:

### AUC_Loss.py
- **Description**: Allows visualization of the loss and AUC (Area Under the Curve) for both training and validation phases of each model training session.
- **Input**: history.pickle file
- **Output**: png with two plot (AUC and Loss)

### Embedding.py
- **Description**: Script to perform an embedding of the weights of the last dense layer of the model. Dimensionality reduction is then performed, allowing visualization of the different samples colored by (patient, cell lines, or hormonal status) in different plots (png).
- **Input**: Count matrix (h5 file), trained model, cvs with hormone status information
- **Output**: 3 embedding PNG plots (leiden, patients vs cell lines, hormone status)

### Overlappin.py
- **Description**: Script that calculates the percentage of overlap between chromosome positions (extended to 1344 bp) within a BED file. 
Possibility to visualize the global overlap distribution through a histogram as well as a chromosome-by-chromosome visualization. 
- **Input**: BED file
- **Output**: overlapping percentage, two png plots




