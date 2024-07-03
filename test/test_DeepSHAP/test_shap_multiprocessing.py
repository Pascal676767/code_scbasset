import os
import pytest
import numpy as np
import sys 
sys.path.append('../../source/DeepSHAP/')
from  shap_functions import *


DATA_FILE = 'test_seqs.h5'



def test_open_file():
    result = open_file(DATA_FILE)
    assert result.shape == (1172, 1344) 


def test_one_hot_encode():
    mock_sequences = np.array([[1, 0, 0, 3, 1, 0, 2, 0, 1]])
    result = one_hot_encode(mock_sequences)
    expected_result = np.array([[[0, 1, 0, 0],
                                 [1, 0, 0, 0],
                                 [1, 0, 0, 0],
                                 [0, 0, 0, 1],
                                 [0, 1, 0, 0],
                                 [1, 0, 0, 0],
                                 [0, 0, 1, 0],
                                 [1, 0, 0, 0],
                                 [0, 1, 0, 0]]], dtype=np.float32)
    assert np.array_equal(result, expected_result), f"Expected {expected_result}, but got {result}"

def test_shuffle_several_times():
    mock_sequences_one_hot = np.array([[[0, 1, 0, 0],
                                 [1, 0, 0, 0],
                                 [1, 0, 0, 0],
                                 [0, 0, 0, 1],
                                 [0, 1, 0, 0],
                                 [1, 0, 0, 0],
                                 [0, 0, 1, 0],
                                 [1, 0, 0, 0],
                                 [0, 1, 0, 0]]], dtype=np.float32)
    result = shuffle_several_times(mock_sequences_one_hot)
    assert result.shape == (100, 9, 4)



if __name__ == '__main__':
    pytest.main()
