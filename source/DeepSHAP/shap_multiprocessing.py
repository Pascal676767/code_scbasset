import h5py
import numpy as np
import shap
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

###############################################################################
# removes tensorflow and shap warning messages, leaving only potential errors #
###############################################################################

logging.getLogger('tensorflow').setLevel(logging.ERROR)
logging.getLogger('shap').setLevel(logging.ERROR)
warnings.filterwarnings('ignore', category=UserWarning)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # '0' = all info, '1' = filter warnings, '2' = filter info, '3' = filter everything except errors

########
# Note #
########

"""
DeepShap is memory-intensive, so the number of cores used for parallelization is greatly limited, 
and is also impacted by the number of reference sequences used.
On the dataset I used, parallelization with 5 cores and 100 shuffles per input sequence requires 
around 200 gigabytes of memory in total. Multiply this by 2 if, for example, you want to use 10 cores. 
The number of sequences to be analyzed doesn't seem to have any impact on memory requirements,
but it does logically have an impact on computation time.
"""




#############
# fonctions #
#############


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
    """
    Input: file path in h5 format
    Output: data matrix contained in file table X 
    """

    with h5py.File(file_path, 'r') as f:
        data = f['X'][:]
        return data

def shap_values(seqs_to_explain, model_path):
    """
    Input: DNA sequence one hot encoded and model path
    Output: DNA sequence with shapley values and DNA sequence one hot encoded
    """
    model = load_model(model_path)

    # the reference to explain the sequences
    dinuc_shuff_explainer = shap.DeepExplainer(model, shuffle_several_times)

    seqs_to_explain = np.expand_dims(seqs_to_explain, axis=0)

    # calculate the shapley value each nucleotide for the input sequence, compare to the reference
    raw_shap_explanations = dinuc_shuff_explainer.shap_values(seqs_to_explain, check_additivity=False)

    # averages all shapley values at a given position 
    dinuc_shuff_explanations = (np.sum(raw_shap_explanations, axis=-1)[:, :] * seqs_to_explain)

    return dinuc_shuff_explanations, seqs_to_explain

def save_incremental_npz(shapley_values, seqs_to_explain, shapley_filename, seqs_filename):
    """
    Input: an encoded one-hot sequence and the shapley values corresponding to this sequence 
    Output: 2 npz files. One for the sequence and one for its shapley values.
    """
    # Save shapley values incrementally
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
    
    np.savez(seqs_filename, combined_seqs)

def task(args):
    """
    Function used for parallelization which uses the shap_value function 
    to calculate the shapley values of the input sequence. 
    It then formats them in the correct dimension so that they can be saved 
    with the save_incremental_npz function. 
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

def shap_sequence_analysis(data_path, model_path, process,  output_folder, data_size):
    """
    Input : 
    - Path to the dataset containing the sequences to be analyzed
    - Path to trained model
    - Number of tasks to be run in parallel (see note)
    - Output folder containing npz files
    - Number of sequences to be analyzed 
    """

    # Data
    data = open_file(data_path)
    data_one_hot = one_hot_encode(data[:data_size])

    with multiprocessing.Pool(process) as pool:
        i = 1
        for shapley_values, seqs_to_explain in tqdm(pool.imap_unordered(task, [(seq, model_path, output_folder) for seq in data_one_hot], chunksize=1), total=len(data_one_hot), leave=True):
            print('Sequence:', i)
            i += 1

            gc.collect()

if __name__ == "__main__":
    multiprocessing.set_start_method('forkserver')
    shap_sequence_analysis(data_path=sys.argv[1], model_path=sys.argv[2], process=int(sys.argv[3]), data_size=int(sys.argv[4]), output_folder=sys.argv[5])
