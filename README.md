# code_scbasset

# Architecture


# Source Folder



## Preprocessing Folder

This folder contains a number of python scripts for pre-processing data in order to train the model. Below you will find descriptions of each notebook:

### CAGE_preprocessing.py
- **Description**: Converts an original CSV file (row: DNA position, column: sample), containing enhancer and promoter data, into a BED file listing the promoter or enhancer regions Also generates an H5 file containing the count matrix and an other H5 file with the coutn matrx filtered (rate to defined).
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



