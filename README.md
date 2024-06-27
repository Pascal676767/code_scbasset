# QUICK START

### Preprocessing

The starting csv file must be in this form (row: DNA position with start, end and strand, column: the various samples)

| Position | Sample 1 | Sample 2 | ........ | Sample n |
|----------|----------|----------|----------|----------|
| chr:xxxx-xxxx:+ | activity | activity | ........ | activity |
| chr:xxxx-xxxx:+ | activity | activity | ........ | activity |
| chr:xxxx-xxxx:- | activity | activity | ........ | activity |


Then use the following command:

```
python CAGE_preprocessing.py [CSV.file], [filter rate] --option
```

You should now have 3 different files: cage.bed, count_matrix.h5, count_matrix_filtered.h5

You now need to download and install the modified version of scBasset available [here](https://github.com/Pascal676767/scBasset).
Refer to part 3. of the Readme of the modified version of scBasset to generate training data from the count_matrix_filtered.h5 files.

### Training

Refer to part 4. of the Readme of the modified version of scBasset to train the model.

### Analysis

Once the model has been trained, use the following command on the history.pickle file:
```
python AUC_Loss.py [history.pickle] --option
```
You should now have an AUC_Loss.png file that lets you view the AUC and Loss on the training and validation data.

To access the embedding view, use the following command:
```
python embedding.py [count_matrix_filtered.h5], [trained model], [csv_hormone]
```
You should now have 3 png files (leiden, patients vs cell lines, hormone status)

### DeepSHAP

To extract what the model has learned from the sequence, we use the DeepSHAP tool. Note that it requires a lot of memory. See note on the script 'shap_multiprocessing.py'.
Use the following command
```
python shap_multiprocessing.py [training, validation or test sequence], [trained model], [core number for parallelization (see note)], [start seq number], [end seq number].
```
You should now have two npz files: seqs_to_explain.npz and shapley_values.npz

### tfModisco-lite

You can now use [tfModisco-lite](https://github.com/jmschrei/tfmodisco-lite) with both npz files as input.



# CODE

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




