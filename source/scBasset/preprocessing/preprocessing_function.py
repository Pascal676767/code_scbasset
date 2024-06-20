import csv
import pandas as pd
from itertools import islice
import h5py
import numpy as np
from scipy import sparse
import re
import scanpy as sc

def adjust_single_coordinates(position):
    """
    If the end is missing in the coordinate, add the start as the end.
    
    Takes the column with the chromosomal position as input.
    Adds a column to the data frame with the correction of missing values.
    
    """
    full_format_regex = r'(chr[\w]+):(\d+)-(\d+)(:[\+\-])?'

    # If the position matches the description.
    match = re.match(full_format_regex, position)
    if match:
        chromosome, start, end, strand = match.groups()
        # If the strand is specified, return the original position.
        if strand is not None:
            return position
        else:
            # Otherwise, assume it's "+".
            return f"{chromosome}:{start}-{end}:+"
    
    # Else, add the start as the end.
    chromosome, start, strand = re.match(r'(chr[\w]+):(\d+):?(\+|-)?', position).groups()
    # If no strand is detected, assume it's "+".
    if strand is None:
        strand = '+'
    return f"{chromosome}:{start}-{start}:{strand}"



def adjust_coordinates(data):
    """
    Input: data frame
    Output: data frame with adjusted coordinates
    """
    df_adjust_coordinates = pd.DataFrame(data)
    # Insert the new column.
    df_adjust_coordinates.insert(1, 'adjusted_coordinate', data['Region'].apply(adjust_single_coordinates))
    return df_adjust_coordinates



def create_bed_file(df_peak_position, output_file='peaks.bed'):
    """
    Creates a BED file from peak coordinates in a pandas DataFrame.

    Args:
    - df_peak_position (pd.Series): A pandas Series containing peak coordinates.
    - output_file (str): The output BED file name. Default is 'peaks.bed'.
    """
    # Extract the information of chromosome, start, end.
    bed_data = df_peak_position.str.extract(r'(chr[\w]+):(\d+)-(\d+):?(\+|-)?', expand=True)

    # Write into a BED file.
    bed_data.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"Successfully created BED file: {output_file}")

def number_features_info(ad):
    print(f'Number of features: {ad.shape[1]} \nNumber of patient: {ad.shape[0]}')


def extract_info(data, df_peak_position):
    # Extract the barcodes (patient IDs), features (peaks), and count matrix from the CSV file.
    # Barcodes
    barcodes = data.columns[1:]
    barcodes_list = list(barcodes.astype(str))

    # Features
    features = df_peak_position
    features_list = list(features.astype(str))
    num_features = len(features_list)
    feature_types_data = np.array(['Peaks'] * len(features_list), dtype='S20')
    genome_data = np.array(['NA'] * len(features_list), dtype='S20')
    
    # Count matrix
    data_matrix = data.iloc[:, 1:].values
    data_matrix_transposed = np.transpose(data_matrix)  # Barcodes (patients) in rows and peaks in columns
    number_features_info(data_matrix_transposed)

    return barcodes_list, data_matrix_transposed, features_list, feature_types_data, genome_data


def create_count_matrix_h5(barcodes_list, data_matrix_transposed, features_list, feature_types_data, genome_data, output_file):
    """
    Creates an HDF5 file with specified groups and datasets for count matrix.

    Args:
    - barcodes_list (list): List of barcodes (patient)
    - data_matrix_transposed: Transposed count matrix data with informations of the activities of the enhancer/promoter
    - features_list (list): List of feature (peaks position)
    - feature_types_data: Type of the features
    - genome_data : Data of genome.
    - output_file (str): Output HDF5 file name.
    """
    with h5py.File(output_file, 'w') as f:
        # Create the /matrix group
        grp_matrix = f.create_group('matrix')
        
        # Write the barcodes into the /matrix/barcodes dataset
        grp_matrix.create_dataset('barcodes', data=barcodes_list)
        
        # Convert the count matrix to CSR format
        data_csr = sparse.csr_matrix(data_matrix_transposed)

        # Write the count matrix CSR into the /matrix/data dataset
        grp_matrix.create_dataset('data', data=data_csr.data)
        
        # Invert the dimensions of the matrix
        data_matrix_shape = data_csr.shape
        data_matrix_shape_inverse = (data_matrix_shape[1], data_matrix_shape[0])
        
        # Write the inverted shape into the /matrix/shape dataset
        grp_matrix.create_dataset('shape', data=data_matrix_shape_inverse)
        
        # Write the indices into the /matrix/indices dataset
        grp_matrix.create_dataset('indices', data=data_csr.indices)
        
        # Write the pointers into the /matrix/indptr dataset
        grp_matrix.create_dataset('indptr', data=data_csr.indptr)
        
        # Create the /matrix/features group
        grp_features = grp_matrix.create_group('features')
        
        # Write information about the features
        grp_features.create_dataset('_all_tag_keys', data=features_list)
        grp_features.create_dataset('feature_type', data=feature_types_data)
        grp_features.create_dataset('genome', data=genome_data)
        grp_features.create_dataset('id', data=features_list)
        grp_features.create_dataset('interval', data=features_list)
        grp_features.create_dataset('name', data=features_list)

        print(f"Successfully created {output_file} file")


def Visualize_h5_file(file_name, group):
    """
    Visualize the contents of an HDF5 file.

    Parameters:
    - file_name (str): The name of the HDF5 file.
    - group (str): The name of the group within the HDF5 file.

    """
    with h5py.File(file_name, 'r') as f:
        # Retrieve the keys of the HDF5 file
        keys = list(f.keys())
        print(f"Keys available in the HDF5 file: {keys}\n")

        # Access the specified group
        group_data = f[group]

        # Iterate over each member of the group
        for member in group_data.keys():
            # Check if the member is a group or a dataset
            if isinstance(group_data[member], h5py.Group):
                # If it's a group, display its keys
                print(f"Keys of subgroup '{member}':", list(group_data[member].keys()))
            else:
                # Otherwise, display the dataset's value
                dimensions = group_data[member].shape
                print(f"Dimension: {dimensions}")
                print(f"Value of dataset '{member}':", group_data[member][:])


def shuffle(csv_file_path):
    df = pd.read_csv(csv_file_path)
    # Exclude column names
    data_matrix = df.iloc[:, 1:].values
    # Random shuffle
    np.random.shuffle(data_matrix)
    df.iloc[:, 1:] = data_matrix
    return df




def filtering(count_matrix_file, bed_file, filter_rate):

    # read count matrix and bed file (peaks)
    peak = pd.read_csv(bed_file, sep='\t', names=['chr','start','end', 'strand'])
    ad = sc.read_10x_h5(count_matrix_file, gex_only=False)
    print('Befor filtering')
    number_features_info(ad)

    ######## filtering of peaks #########
    ad_cage = ad[:, ad.var['feature_types']=='Peaks']
    ad_cage.var['chr'] = peak['chr'].values
    ad_cage.var['start'] = peak['start'].values
    ad_cage.var['end'] = peak['end'].values
    ad_cage.var['strand'] = peak['strand'].values

    sc.pp.filter_cells(ad_cage, min_genes=0)
    sc.pp.filter_genes(ad_cage, min_cells=0)

    # a peak need to be accessible in _% cells
    thres = int(ad.shape[0]*filter_rate)
    ad_cage = ad_cage[:, ad_cage.var['n_cells']>thres]
    print('After peak filtering')
    number_features_info(ad_cage)

    ######## filtering of chromosome ########
    chrs = ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
    ad_cage = ad_cage[:, ad_cage.var['chr'].isin(chrs)]
    print('After chr filtering')
    number_features_info(ad_cage)
    
    return ad_cage


def concat_csv_file(file1, file2):
    """
    Input: two csv file
    Output: merge csv file
    """
    data_file1 = pd.read_csv(file1, sep = ',', header = 0)
    data_file2 = pd.read_csv(file2, sep = ',', header = 0)
    
    if (data_file1.columns == data_file2.columns).all():
        data_concat = pd.concat([data_file1, data_file2], ignore_index=True)
        return data_concat
        








