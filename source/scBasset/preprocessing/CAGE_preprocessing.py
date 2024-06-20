import csv
import pandas as pd
from itertools import islice
import h5py
import numpy as np
from scipy import sparse
import re
import scanpy as sc
import sys
from preprocessing_function import *


import csv
import pandas as pd
from itertools import islice
import h5py
import numpy as np
from scipy import sparse
import re
import scanpy as sc
import sys
from preprocessing_function import *


def CAGE_preprocessing(csv_file, output_bed_file, output_count_matrix, output_count_matrix_filtered, filter_rate):
    """
    Input: CSV file (row: DNA position, column: sample)
    Output: Bed file with the DNA position, count matrix (h5 file) and count matrix filtered (h5 file)
    """

    # Load the data
    data = open_verif_csv(csv_file)

    # Completing the chromosomal coordinates (make sure every DNA position has a start and an end)
    df_adjust_coordinates = adjust_coordinates(data)
    df_peak_position = df_adjust_coordinates['adjusted_coordinate']

    # Create the bed file
    create_bed_file(df_peak_position, output_bed_file)

    # Extract the necessary information to create the h5 file
    barcodes_list, data_matrix_transposed, features_list, feature_types_data, genome_data = extract_info(data, df_peak_position)

    # Create the h5 file
    create_count_matrix_h5(barcodes_list, data_matrix_transposed, features_list, feature_types_data, genome_data, output_count_matrix)

    # Create the filtered h5 file
    ad_cage = filtering(output_count_matrix, output_bed_file, filter_rate)
    ad_cage.write(output_count_matrix_filtered)

def print_help():
    help_text = """
    Usage: python script.py <csv_file> <output_bed_file> <output_count_matrix> <output_count_matrix_filtered> <filter_rate>
    
    Arguments:
    csv_file                  Path to the input CSV file.
    output_bed_file           Path to the output BED file.
    output_count_matrix       Path to the output count matrix H5 file.
    output_count_matrix_filtered Path to the output filtered count matrix H5 file.
    filter_rate               Filtering rate for the count matrix.
    
    Options:
    --help                    Show this help message and exit.
    """
    print(help_text)

if __name__ == "__main__":
    if '--help' in sys.argv or len(sys.argv) != 6:
        print_help()
    else:
        csv_file = sys.argv[1]
        output_bed_file = sys.argv[2]
        output_count_matrix = sys.argv[3]
        output_count_matrix_filtered = sys.argv[4]
        filter_rate = float(sys.argv[5])  # Convert the filter rate to a float
        CAGE_preprocessing(csv_file, output_bed_file, output_count_matrix, output_count_matrix_filtered, filter_rate)
