import pandas as pd
import numpy as np
import sys
from preprocessing_function import shuffle


def controle_shuffle(csv_file_path, output_csv_file_path):
    """
    Input: CSV file (row: DNA position, column: sample)
    Output: Shuffled CSV data file
    """
    df_shuffle = shuffle(csv_file_path)
    df_shuffle.to_csv(output_csv_file_path, index=False)


def print_help():
    """
    Prints the help message for the script usage.
    """
    help_text = """
    Usage: python script.py <csv_file> <output_csv_file>
    
    Arguments:
    csv_file           Path to the input CSV file.
    output_csv_file    Path to the output shuffled CSV file.
    
    Options:
    --help             Show this help message and exit.
    """
    print(help_text)



if __name__ == "__main__":
    if '--help' in sys.argv or len(sys.argv) != 3:
        print_help()
    else:
        csv_file_path = sys.argv[1]
        output_csv_file_path = sys.argv[2]
        controle_shuffle(csv_file_path, output_csv_file_path)
