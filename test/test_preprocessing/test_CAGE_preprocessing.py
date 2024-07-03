import sys 
sys.path.append('../source/scBasset/preprocessing/')
import os

from preprocessing_function import *


def test_adjust_coordinates_good():
    test_good_position = 'chr1:1-2:+'
    actual_position = adjust_single_coordinates(test_good_position)
    expected_position = 'chr1:1-2:+'
    
    assert expected_position == actual_position

def test_adjust_coordinates_bad():
    test_bad_position = 'chr1:1:-'
    actual_position = adjust_single_coordinates(test_bad_position)
    expected_position = 'chr1:1-1:-'
    
    assert expected_position == actual_position


def test_extract_info():

    data = pd.DataFrame({
    'Region': ['chr1:586071-586624', 'chr1:827333-827837', 'chr1:899070-899511'],
    's15t1p2sq1': [0, 40, 0],
    's15t2p2sq1': [2, 135, 0],
    's15t3p2sq1': [0, 186, 0],
    })

    df_peak_position = data['Region']


    barcodes_list, data_matrix_transposed, features_list, feature_types_data, genome_data = extract_info(data, df_peak_position)
    
    expected_barcodes = ['s15t1p2sq1', 's15t2p2sq1', 's15t3p2sq1'] 
    expected_data_matrix_transposed = np.array([
        [0, 40, 0],
        [2, 135, 0],
        [0, 186, 0]
    ])
    expected_features_list = ['chr1:586071-586624', 'chr1:827333-827837', 'chr1:899070-899511']
    expected_feature_types_data = np.array(['Peaks'] * 3, dtype='S20')
    expected_genome_data = np.array(['NA'] * 3, dtype='S20') 

    assert barcodes_list == expected_barcodes
    assert np.array_equal(data_matrix_transposed, expected_data_matrix_transposed)
    assert features_list == expected_features_list
    assert np.array_equal(feature_types_data, expected_feature_types_data)
    assert np.array_equal(genome_data, expected_genome_data)
    

def test_create_count_matrix():
    barcodes_list = ['s15t1p2sq1', 's15t2p2sq1', 's15t3p2sq1'] 
    data_matrix_transposed = np.array([
        [0, 40, 0],
        [2, 135, 0],
        [0, 186, 0]
    ])
    features_list = ['chr1:586071-586624', 'chr1:827333-827837', 'chr1:899070-899511']
    feature_types_data = np.array(['Peaks'] * 3, dtype='S20')
    genome_data = np.array(['NA'] * 3, dtype='S20')
    output_file = './test_create_count_matrix.h5'

    create_count_matrix_h5(barcodes_list, data_matrix_transposed, features_list, feature_types_data, genome_data, output_file) 

    with h5py.File(output_file, 'r') as f:
        assert 'matrix' in f
        assert 'features' in f['matrix']

        # Change bytes in str
        hdf5_barcodes = f['matrix/barcodes'][:]
        hdf5_barcodes_str = [barcode.decode('utf-8') for barcode in hdf5_barcodes]

        hdf5_features = f['matrix/features/_all_tag_keys'][:]
        hdf5_features_str = [feature.decode('utf-8') for feature in hdf5_features]

        hdf5_id = f['matrix/features/id'][:]
        hdf5_id_str = [id.decode('utf-8') for id in hdf5_id]

        hdf5_interval = f['matrix/features/interval'][:]
        hdf5_interval_str = [interval.decode('utf-8') for interval in hdf5_interval]

        hdf5_name = f['matrix/features/name'][:]
        hdf5_name_str = [name.decode('utf-8') for name in hdf5_name]


        assert np.array_equal(hdf5_barcodes_str, barcodes_list)
        assert np.array_equal(f['matrix/data'][:], sparse.csr_matrix(data_matrix_transposed).data)
        assert np.array_equal(f['matrix/shape'][:], np.array([3, 3]))
        assert np.array_equal(f['matrix/indices'][:], sparse.csr_matrix(data_matrix_transposed).indices)
        assert np.array_equal(f['matrix/indptr'][:], sparse.csr_matrix(data_matrix_transposed).indptr)
        assert np.array_equal(hdf5_features_str, features_list)
        assert np.array_equal(f['matrix/features/feature_type'][:], feature_types_data)
        assert np.array_equal(f['matrix/features/genome'][:], genome_data)
        assert np.array_equal(hdf5_id_str, features_list)
        assert np.array_equal(hdf5_interval_str, features_list)
        assert np.array_equal(hdf5_name_str, features_list)

    os.remove(output_file)





