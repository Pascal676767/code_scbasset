import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
import h5py


def read_mtx(matrix_file):
    """
    This function reads the sparse matrix file.
    Input: sparse matrix file in mtx format
    Output: number of rows, columns, and the data of the sparse matrix
    """
    with open(matrix_file, 'r') as f:
        # Skip headers
        for _ in range(2):
            next(f)
        data, row, col = [], [], []
        for line in f:
            r, c, d = map(int, line.split())
            row.append(r - 1)  # adjust for zero-based indexing
            col.append(c - 1)  # adjust for zero-based indexing
            data.append(d)
    return row, col, data

def read_names(file):
    with open(file, 'r') as f:
        return [line.strip() for line in f]

def create_h5_from_files(matrix='file.mtx', column='file_column.txt', row='file_row.txt', output_file='final_matrix_transposed.h5'):
    """
    This function creates an h5 file.
    Input: sparse matrix, txt file with column names, txt file with row names
    Output: h5 file
    """
    
    # Call the function to read the sparse matrix
    row_indices, col_indices, data = read_mtx(matrix)

    sparse_matrix = coo_matrix((data, (row_indices, col_indices)))
    csr_matrix_data = sparse_matrix.tocsr()

    # Transpose the CSR matrix because we want the cells in rows and the features (peaks) in columns
    csr_matrix_transposed = csr_matrix_data.transpose().tocsr()

    # Read column and row names
    barcodes_list = np.array(read_names(column), dtype='S')  
    features_list = np.array(read_names(row), dtype='S')  
    
    # Add some information needed to be read by some other scripts
    feature_types_data = np.array(['Peaks'] * len(features_list), dtype='S')
    genome_data = np.array(['GRCh38'] * len(features_list), dtype='S')
    
    with h5py.File(output_file, 'w') as f:
        # Create the /matrix group
        grp_matrix = f.create_group('matrix')
        
        # Write barcodes (column names in this context)
        grp_matrix.create_dataset('barcodes', data=barcodes_list)
        
        # Write CSR components of the transposed count matrix
        grp_matrix.create_dataset('data', data=csr_matrix_transposed.data)
        grp_matrix.create_dataset('indices', data=csr_matrix_transposed.indices)
        grp_matrix.create_dataset('indptr', data=csr_matrix_transposed.indptr)

        # Very weird, but needs to be the inverse of the actual shape to be read by the next script implemented in scbasset.
        grp_matrix.create_dataset('shape', data=csr_matrix_data.shape)
        
        # Create the /matrix/features group
        grp_features = grp_matrix.create_group('features')
        
        # Write feature metadata
        grp_features.create_dataset('_all_tag_keys', data=features_list)
        grp_features.create_dataset('feature_type', data=feature_types_data)
        grp_features.create_dataset('genome', data=genome_data)
        grp_features.create_dataset('id', data=features_list)
        grp_features.create_dataset('interval', data=features_list)
        grp_features.create_dataset('name', data=features_list)

        print(f"The file {output_file} has been successfully created.")

# Calling the function
create_h5_from_files(
    matrix='../../../results_scbasset/data/SCAFE/Malignant_cells_distal_cre.mtx',
    column='../../../results_scbasset/data/SCAFE/colnames_Malignant_cells_distal_cre.txt',
    row='../../../results_scbasset/data/SCAFE/rownames_Malignant_cells_distal_cre.txt',
    output_file='../../../results_scbasset/enhancer_SCAFE/data/count_matrix_enhancer.h5'
)


