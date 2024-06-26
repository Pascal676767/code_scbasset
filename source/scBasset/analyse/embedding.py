from analyse_function import *
import os
import argparse
import anndata


"""
Script to perform an embedding of the weights of the last dense layer of the model. 
Dimensionality reduction is then performed, allowing visualization of the different samples colored by (patient, cell lines, or hormonal status) in different plots (png).
"""

def main_embedding(anndata_file, model_path, hormone_path, results_folder='./results_embedding'):
    """
        Permits to visualize the embedding of the last dense layer weights, colored by patient, cell lines, and hormonal status.

        Parameters:
        - anndata_file: path to the AnnData file containing the samples and the features
        - model_path: path to the trained model
        - hormone_path: path to the CSV file containing the hormonal status of the patients
        - results_folder: folder in which the embeddings will be stored

        Output:
        - 3 embedding PNG plots
    """

    os.makedirs(results_folder, exist_ok=True)

    ad = anndata.read_h5ad(anndata_file)
    model = load_model(ad, model_path)

    embedding(model, results_folder)

    add_patient_info(ad)
    add_hormone_info(ad, hormone_path)

    plotting_embedding(ad, 'leiden', results_folder)
    plotting_embedding(ad, 'Hormone',results_folder)
    plotting_embedding(ad, 'Patient', results_folder)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('anndata_file', type=str, help='Path to the anndata file containing the DNA position and features.')
    parser.add_argument('model_path', type=str, help='Path to the trained model')
    parser.add_argument('hormone_path', type=str, help='Path to the hormone file')
    parser.add_argument('--results_folder', type=str, default='./results_embedding', help='Folder path to save the plots.')

    args = parser.parse_args()

    main_embedding(anndata_file=args.anndata_file, model_path=args.model_path, hormone_path=args.hormone_path, results_folder=args.results_folder)


if __name__ == '__main__':
    main()
