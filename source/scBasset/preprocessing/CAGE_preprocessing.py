import argparse
from preprocessing_function import *
import warnings

# Ignore all warnings
warnings.filterwarnings('ignore')



def CAGE_preprocessing(csv_file, filter_rate, output_bed_file='./cage.bed', output_count_matrix='./count_matrix.h5', output_count_matrix_filtered='./count_matrix_filtered.h5'):
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

def main():
    parser = argparse.ArgumentParser(description='Preprocess CAGE data from a CSV file and generate BED and H5 files.')
    parser.add_argument('csv_file', type=str, help='Path to the input CSV file.')
    parser.add_argument('filter_rate', type=float, help='Filtering rate for the count matrix.')
    parser.add_argument('--output_bed_file', type=str, default='./cage.bed', help='Path to the output BED file. Default is "./cage.bed".')
    parser.add_argument('--output_count_matrix', type=str, default='./count_matrix.h5', help='Path to the output count matrix H5 file. Default is "./count_matrix.h5".')
    parser.add_argument('--output_count_matrix_filtered', type=str, default='./count_matrix_filtered.h5', help='Path to the output filtered count matrix H5 file. Default is "./count_matrix_filtered.h5".')


    args = parser.parse_args()

    CAGE_preprocessing(
        csv_file=args.csv_file,
        filter_rate=args.filter_rate,
        output_bed_file=args.output_bed_file,
        output_count_matrix=args.output_count_matrix,
        output_count_matrix_filtered=args.output_count_matrix_filtered
    )

if __name__ == "__main__":
    main()
