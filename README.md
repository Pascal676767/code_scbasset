# code_scbasset

# Source Folder

The source folder contains two main folder and several bash script

## Preprocessing Folder

This folder contains a number of jupyter notebooks for pre-processing data in order to train the model. Below you will find descriptions of each notebook:

### CAGE_preprocessing.ipynb
- **Description**: Converts an original CSV file (containing enhancer and promoter data) into a BED file listing the peak regions (promoter and enhancer). Also generates an H5 file containing the count matrix.
- **Input**: CSV file (enhancer and promoter).
- **Output**:
  - BED file of peaks.
  - H5 file (count matrix).

### Concat_enhancer_promoteur.ipynb
- **Description**: Concatenates two CSV files, typically those of promoters and enhancers.
- **Input**: Two CSV files (promoter and enhancer).
- **Output**: A combined CSV file.

### Controle_shuffle.ipynb
- **Description**: Randomly shuffles the values in the original CSV file to create a control.
- **Input**: Original CSV file.
- **Output**: Shuffled CSV file.

### SCAFE_preprocessing.py
- **Description**: From a sparse matrix and two text files containing column and row information, creates an H5 file containing the count matrix.
- **Input**:
  - Sparse matrix.
  - Text files for columns and rows.
- **Output**: H5 file (count matrix).

### Testing-random_stuff.ipynb
- **Description**: A no-op script intended for various tests.
- **Function**: None (testing purposes only).

### make_anndata_CAGE.ipynb
- **Description**: Generates the final file necessary for creating training files. Allows filtering of the number of enhancers, promoters, etc.
- **Function**: Data preparation for training.


## Analysis Folder

The analysis directory includes several Jupyter Notebooks designed for visualizing various aspects of model training and evaluation. Below you will find descriptions of each notebook:

### AUC_Loss.ipynb
- **Description**: Allows visualization of the loss and AUC (Area Under the Curve) for both training and validation phases of each model training session.
- **Purpose**: To evaluate model performance and overfitting by comparing training and validation metrics.

### Embedding_Pearson.ipynb
- **Description**: Facilitates the visualization of the embedding of patient vs. cell line for hormones.
- **Purpose**: To understand the relationships and similarities between patient and cell line hormone data through embedding visualizations.

### Probability_visualisation.ipynb
- **Description**: Enables visualization of the class proportions (0 and 1) and the probability distribution of train and validation values for the loss.
- **Purpose**: To assess the model's ability to classify and predict by examining the distribution of probabilities and class proportions.

### Score_motif.ipynb
- **Description**: Allows for the visualization of the activity of a particular motif.
- **Purpose**: To analyze and interpret the significance and activity of specific motifs in the context of the model's predictions.


