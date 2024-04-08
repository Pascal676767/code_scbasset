import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import h5py

# Define functions to read each type of file
def read_mtx(matrix_file):
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
    # Read matrix data
    row_indices, col_indices, data = read_mtx(matrix)
    sparse_matrix = coo_matrix((data, (row_indices, col_indices)))
    csr_matrix_data = sparse_matrix.tocsr()
    
    # Transpose the CSR matrix
    csr_matrix_transposed = csr_matrix_data.transpose()
    
    # Read column and row names
    barcodes_list = np.array(read_names(column), dtype='S')  
    features_list = np.array(read_names(row), dtype='S')  
    
 
    feature_types_data = np.array(['gene_expression'] * len(features_list), dtype='S')
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
        grp_matrix.create_dataset('shape', data=csr_matrix_transposed.shape)
        
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

# Example of calling the function
create_h5_from_files(
    matrix='../../../results_scbasset/data/SCAFE/Malignant_cells_distal_cre.mtx',
    column='../../../results_scbasset/data/SCAFE/colnames_Malignant_cells_distal_cre.txt',
    row='../../../results_scbasset/data/SCAFE/rownames_Malignant_cells_distal_cre.txt',
    output_file='../../../results_scbasset/enhancer_SCAFE/data/count_matrix_enhancer.h5'
)


