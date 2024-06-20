import gzip
import h5py
import numpy as np
import shap
from deeplift.visualization import viz_sequence
from keras.models import load_model
shap.explainers._deep.deep_tf.op_handlers["AddV2"] = shap.explainers._deep.deep_tf.passthrough
from deeplift.dinuc_shuffle import dinuc_shuffle
import sys
import multiprocessing
from tqdm import tqdm
import gc

import logging
import os
import warnings

logging.getLogger('tensorflow').setLevel(logging.ERROR)
logging.getLogger('shap').setLevel(logging.ERROR)

warnings.filterwarnings('ignore', category=UserWarning)


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # '0' = toutes les infos, '1' = filtrer les warnings, '2' = filtrer les infos, '3' = filtrer tout sauf les erreurs


def one_hot_encode(sequences):
    """
    Input: DNA sequence
    Output: DNA sequence onehot encoded
    """
    one_hot = np.zeros((sequences.shape[0], sequences.shape[1], 4), dtype=np.float32)

    for i in range(4):
        one_hot[:, :, i] = (sequences == i)

    return one_hot

def shuffle_several_times(s):
    """
    Input: DNA sequence one hot encoded
    Output: DNA sequence X time shuffeld
    """
    s = np.squeeze(s)
    return dinuc_shuffle(s, num_shufs=100)

def open_file(file_path):
    with h5py.File(file_path, 'r') as f:
        x_test = f['X'][:]
        return x_test

def shap_values(seqs_to_explain, model_path):
    """
    Input: DNA sequence one hot encoded
    Output: DNA sequence with shapley values and DNA sequence one hot encoded
    """
    model = load_model(model_path)

    # the reference to explain the sequences
    dinuc_shuff_explainer = shap.DeepExplainer(model, shuffle_several_times)

    seqs_to_explain = np.expand_dims(seqs_to_explain, axis=0)

    # calculate the shapley value each nucleotide for the input sequence, compare to the reference
    raw_shap_explanations = dinuc_shuff_explainer.shap_values(seqs_to_explain, check_additivity=False)

    dinuc_shuff_explanations = (np.sum(raw_shap_explanations, axis=-1)[:, :] * seqs_to_explain)

    return dinuc_shuff_explanations, seqs_to_explain

def save_incremental_npz(shapley_values, seqs_to_explain, shapley_filename, seqs_filename):
    """
    Save the data incrementally into two NPZ files.
    """
    # Save shapley values incrementally
    print('Writing')
    try:
        existing_shapley = np.load(shapley_filename)['arr_0']
        combined_shapley = np.concatenate((existing_shapley, shapley_values), axis=0)
    except FileNotFoundError:
        combined_shapley = shapley_values

    np.savez(shapley_filename, combined_shapley)

    # Save seqs to explain incrementally
    try:
        existing_seqs = np.load(seqs_filename)['arr_0']
        combined_seqs = np.concatenate((existing_seqs, seqs_to_explain), axis=0)
    except FileNotFoundError:
        combined_seqs = seqs_to_explain
    
    print('Writing done')


    np.savez(seqs_filename, combined_seqs)

def task(args):
    """
    Input: Seqs to explain
    Output: Seqs to explain with shapley values
    """
    i, model, output_folder = args
    shapley_values, seqs_to_explain = shap_values(i, model)

    shapley_values = np.squeeze(shapley_values)
    seqs_to_explain = np.squeeze(seqs_to_explain)

    shapley_values = np.expand_dims(shapley_values, axis=0)
    seqs_to_explain = np.expand_dims(seqs_to_explain, axis=0)

    shapley_filename = f"{output_folder}/shapley_values.npz"
    seqs_filename = f"{output_folder}/seqs_to_explain.npz"

    save_incremental_npz(shapley_values, seqs_to_explain, shapley_filename, seqs_filename)

    return shapley_values, seqs_to_explain

def data_modisco(data_path, model_path, process,  output_folder, data_size):
    # Data
    data = open_file(data_path)
    data_one_hot = one_hot_encode(data[:data_size])

    num_processes = multiprocessing.cpu_count()
    print('nombre de processus disponible ', num_processes)

    # shapley_filename = f"{output_folder}/shapley_values.npz"
    # seqs_filename = f"{output_folder}/seqs_to_explain.npz"

    with multiprocessing.Pool(process) as pool:
        i = 1
        for shapley_values, seqs_to_explain in tqdm(pool.imap_unordered(task, [(seq, model_path, output_folder) for seq in data_one_hot], chunksize=1), total=len(data_one_hot), leave=True):
            print('Sequence:', i)
            # Save results incrementally in the NPZ files
            # save_incremental_npz(shapley_values, seqs_to_explain, shapley_filename, seqs_filename)
            i += 1

            gc.collect()

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.set_start_method('forkserver')
    data_modisco(data_path=sys.argv[1], model_path=sys.argv[2], process=int(sys.argv[3]), data_size=int(sys.argv[4]), output_folder=sys.argv[5])
