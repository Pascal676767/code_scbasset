import anndata
import tensorflow as tf
import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import math
import pickle
import seaborn as sns
import scipy
import sys
import scanpy as sc
import pandas as pd
from scbasset.utils import *





def load_model(ad, trained_model):
    """
    Input: anndata file and model weight
    """
    # load model
    model = make_model(32, ad.shape[0], show_summary=False)
    model.load_weights(trained_model)
    return model

def intercept_depth(ad, model, type_data):
    intercept = get_intercept(model) # get_intercept function
    sc.pp.filter_cells(ad, min_counts=0)

    f, ax = plt.subplots(figsize=(4,4))
    r = scipy.stats.pearsonr(intercept, np.log10(ad.obs['n_genes']))[0]
    sns.scatterplot(x=intercept, y=np.log10(ad.obs['n_genes']), ax=ax)
    ax.set_xlabel('intercept')
    ax.set_ylabel('log10(n_peaks)')
    ax.set_title(f'{type_data}_Pearson R: %.3f'%r)

def embedding(model, type_data, results_path):
    proj = get_cell_embedding(model) # get_cell_embedding function
    pd.DataFrame(proj).to_csv(f'{results_path}/projection_cage_{type_data}.csv')

def plotting_embedding(ad_file, color_by, type_data, results_path, option=False, title=None):
    """
    Input: ad_file
    Output: embedding of the 'color_by'
    """
    f, ax = plt.subplots(figsize=(4, 4))
    
    if option:
        projection_data = pd.read_csv(f'{results_path}/projection_cage_{type_data}.csv', index_col=0, nrows=72)
        ad_file.obsm['projection'] = projection_data.values
    else:
        ad_file.obsm['projection'] = pd.read_csv(f'{results_path}/projection_cage_{type_data}.csv', index_col=0).values
    
    sc.pp.neighbors(ad_file, use_rep='projection')
    sc.tl.umap(ad_file)
    sc.tl.leiden(ad_file)
    
    # Use the specified title if provided, otherwise use the default title
    if title is None:
        title = f'{type_data} embedding'
    
    sc.pl.umap(ad_file, color=color_by, ax=ax, size=350, title=title)


def add_patient_info(ad):
    """
    Input: anndata file
    Output: anndata file with patient and cell line information
    """
    ad.obs['Patient'] = np.where(ad.obs.index.str.startswith('s'), 'patient', 'cell_line')
    ad.obs['Patient'] = ad.obs['Patient'].astype('category')

def hormone_classe(row):
    """
    Input: row of hormone csv
    Output: return the hormone class
    """

    #Luminal A
    if row['Oestrogen'] == 1 and row['progesteron'] == 1 and row['her2'] == 0:
        return 'Luminal A'
    # Luminal B(Her2+)
    elif row['Oestrogen'] == 1 and row['progesteron'] == 1 and row['her2'] == 1:
        return 'Luminal B (HER2+)'
    elif row['Oestrogen'] == 1 and row['progesteron'] == 0 and row['her2'] == 1:
        return 'Luminal B (HER2+)'
    # Luminal B(Her2-)
    elif row['Oestrogen'] == 1 and row['progesteron'] == 0 and row['her2'] == 0:
        return 'Luminal B (HER2-)'
    elif row['Oestrogen'] == 0 and row['progesteron'] == 1 and row['her2'] == 0:
        return 'Luminal B (HER2-)'

    # Her2-enriched
    elif row['Oestrogen'] == 0 and row['progesteron'] == 0 and row['her2'] == 1:
        return 'Her2-enriched'
    # Basal-like
    elif row['Oestrogen'] == 0 and row['progesteron'] == 0 and row['her2'] == 0:
        return 'Basal-like'

def add_hormone_info(ad, hormone_path):
    """
    Input: anndata file
    Output: anndata file with hormone information 
    """
    hormone = pd.read_csv(hormone_path, sep = ',', header = 0)
    hormone['Classe'] = hormone.apply(hormone_classe, axis=1)
    sample_class_map = dict(zip(hormone['Sample_id'], hormone['Classe']))
    ad.obs['Hormone'] = [sample_class_map.get(sample_id, None) for sample_id in ad.obs.index]
