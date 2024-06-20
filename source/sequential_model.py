from deeplift.conversion import kerasapi_conversion as kc
from scbasset.utils import *
from scbasset.basenji_utils import *
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import custom_object_scope
import json
from tensorflow.keras.regularizers import l2

def make_sequential_model(
    bottleneck_size,
    n_cells,
    seq_len=1344,
    show_summary=True,
):

    input_config = {
    "batch_input_shape": [None, seq_len, 4],
    "dtype": "float32",
    "sparse": False,
    "ragged": False,
    "name": "sequence"
    }


    model = tf.keras.models.Sequential()
    #0
    model.add(tf.keras.layers.InputLayer(**input_config))
    # model.add(tf.keras.Input(input_shape=(1344, 4), name='sequence'))
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=288, kernel_size=17,strides=1, padding= 'same', name='conv1d', use_bias=False, kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization'
    ))
    model.add(tf.keras.layers.MaxPooling1D(
        pool_size=3, 
        strides=3, 
        padding="same",
        data_format="channels_last",
        name='max_pooling1d',
    ))
    #1
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=288, kernel_size=5,
                                strides=1, use_bias=False, padding='same', name="conv1d_1", kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_1'
    ))
    model.add(tf.keras.layers.MaxPooling1D(
        pool_size=2, 
        strides=2, 
        padding="same",
        data_format="channels_last",
        name='max_pooling1d_1',
    ))
    #2
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=323, kernel_size=5,
                                strides=1, use_bias=False, padding='same', name="conv1d_2", kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_2'
    ))
    model.add(tf.keras.layers.MaxPooling1D(
        pool_size=2, 
        strides=2, 
        padding="same",
        data_format="channels_last",
        name='max_pooling1d_2',
    ))
    #3
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=363, kernel_size=5,
                                strides=1, use_bias=False,  padding='same',name="conv1d_3", kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_3'
    ))
    model.add(tf.keras.layers.MaxPooling1D(
        pool_size=2, 
        strides=2, 
        padding="same",
        data_format="channels_last",
        name='max_pooling1d_3',
    ))
    #4
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=407, kernel_size=5,
                                strides=1, use_bias=False, padding='same', name="conv1d_4", kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_4'
    ))
    model.add(tf.keras.layers.MaxPooling1D(
        pool_size=2, 
        strides=2, 
        padding="same",
        data_format="channels_last",
        name='max_pooling1d_4',
    ))
    #5
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=456, kernel_size=5,
                                strides=1, use_bias=False, padding='same', name="conv1d_5", kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_5'
    ))
    model.add(tf.keras.layers.MaxPooling1D(
        pool_size=2, 
        strides=2, 
        padding="same",
        data_format="channels_last",
        name='max_pooling1d_5',
    ))
    #6
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=512, kernel_size=5,
                                strides=1, use_bias=False, padding='same', name="conv1d_6", kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_6'
    ))
    model.add(tf.keras.layers.MaxPooling1D(
        pool_size=2,
        strides=2, 
        padding="same",
        data_format="channels_last",
        name='max_pooling1d_6',
    ))
    #7
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Conv1D(filters=256, kernel_size=1,
                                strides=1, use_bias=False, padding='same',  name="conv1d_7", kernel_regularizer=l2(0)))
    model.add(tf.keras.layers.BatchNormalization(axis=2,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_7'
    ))
    #8
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Flatten(name='flatten'))
    model.add(tf.keras.layers.Dense(units=bottleneck_size,use_bias=False, activation= "linear", name='dense'))
    model.add(tf.keras.layers.BatchNormalization(axis=1,
        momentum=0.9,
        epsilon=0.001,
        center=True,
        scale=True,
        beta_initializer="zeros",
        gamma_initializer="ones",
        moving_mean_initializer="zeros",
        moving_variance_initializer="ones",
        beta_regularizer=None,
        gamma_regularizer=None,
        beta_constraint=None,
        gamma_constraint=None,
        synchronized=False,
        name= 'batch_normalization_8'
    ))
    model.add(tf.keras.layers.Dropout(0.5))
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.Dense(units =n_cells, activation='sigmoid', use_bias=True, name='dense_1', kernel_regularizer=l2(0.1)))
    model.add(tf.keras.layers.Flatten(data_format='channels_last', name='flatten_1'))
 
    # model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
    # model.build((None, 1344, 4))

    if show_summary:
        model.summary()
    return model
