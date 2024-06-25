import argparse
import warnings
import logging
import os
import multiprocessing

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



if __name__ == "__main__":
    multiprocessing.set_start_method('forkserver')

    # Argument parsing
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('data_path', type=str, help='Path to the dataset')
    parser.add_argument('model_path', type=str, help='Path to trained model')
    parser.add_argument('process', type=int, help='Number of parallel processes')
    parser.add_argument('start_idx', type=int, help='Start index for data range')
    parser.add_argument('end_idx', type=int, help='End index for data range')
    parser.add_argument('output_folder', type=str, help='Output folder path')

    args = parser.parse_args()

    # Calling function with parsed arguments
    shap_sequence_analysis(
        data_path=args.data_path,
        model_path=args.model_path,
        process=args.process,
        output_folder=args.output_folder,
        data_range=(args.start_idx, args.end_idx)
    )
