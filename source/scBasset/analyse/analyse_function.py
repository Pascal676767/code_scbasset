import anndata
import tensorflow as tf
import numpy as np
import pickle
import scipy
import scanpy as sc
import pandas as pd
from scbasset.utils import *
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D
import re
tf.get_logger().setLevel('ERROR')



#####################################
############# AUC_Loss ##############
#####################################


def plot_metric(train_metric, val_metric, title, ylabel, subplot_position):
    """
        Plot the training and validation metrics.
        
        Parameters:
        - metric: list of training metric values.
        - val_metric: list of validation metric values.
        - title: title of the plot.
        - ylabel: label for the y-axis.
        - subplot_position: position of the subplot in the figure.
    """
    plt.subplot(2, 1, subplot_position)
    plt.plot(train_metric, label='Train', color='orange')
    plt.plot(val_metric, label='Validation', color='blue')
    plt.ylim(0, 1)
    plt.title(title)
    plt.xlabel('Epoch')
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)


def train_validation_comparison(file_path, title_loss ='Comparison between train and validation Loss' , 
                                title_auc = 'Comparison between train and validation AUC', output_file ='./AUC_Loss.png'):
    """
        Plot the AUC and Loss for both validation and training data

        Parameters:
        - file_path: history.pickle file
        - title_loss: title for the Loss plot
        - title_auc: title for the AUC plot
        - output_file: name of the output PNG file with the plots

        Output:
        - PNG file
    """


    # Load the data from the pickle file
    with open(file_path, 'rb') as f:
        data = pickle.load(f)

    # Extract the data
    loss = data['loss']
    val_loss = data['val_loss']
    auc = data['auc']
    val_auc = data['val_auc']

    # Create a plot with two subplots (one for the loss and the other for the AUC)
    plt.figure(figsize=(5, 8)) 

    # Plot for the AUC
    plot_metric(auc, val_auc, title_auc, 'AUC', 1)

    # Plot for the Loss
    plot_metric(loss, val_loss, title_loss, 'Loss', 2)

    plt.tight_layout()

    plt.savefig(output_file, format='png')




#################################
########### Embedding ###########
#################################

def load_model(ad, trained_model):
    """
        Create the model and load the corresponding weights

        Parameters:
        - ad: h5 file with features and samples
        - trained_model: weights of the trained model

        Output:
        - Model with weights
    """

    model = make_model(32, ad.shape[0], show_summary=False)
    model.load_weights(trained_model)
    return model

# def intercept_depth(ad, model, type_data):
#     intercept = get_intercept(model) # get_intercept function
#     sc.pp.filter_cells(ad, min_counts=0)

#     f, ax = plt.subplots(figsize=(4,4))
#     r = scipy.stats.pearsonr(intercept, np.log10(ad.obs['n_genes']))[0]
#     sns.scatterplot(x=intercept, y=np.log10(ad.obs['n_genes']), ax=ax)
#     ax.set_xlabel('intercept')
#     ax.set_ylabel('log10(n_peaks)')
#     ax.set_title(f'{type_data}_Pearson R: %.3f'%r)

def embedding(model, type_data, results_path):
    proj = get_cell_embedding(model) # get_cell_embedding function
    pd.DataFrame(proj).to_csv(f'{results_path}/projection_cage_{type_data}.csv')

def plotting_embedding(ad_file, color_by, type_data, results_path, option = False):
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
    sc.pl.umap(ad_file, color=color_by, ax=ax, size=350, title = f'{type_data} embedding')

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
    if row['ER+_HER2-'] == 1:
        return 'ER+_HER2-'
    elif row['HER2+'] == 1:
        return 'HER2+'
    elif row['TNBC'] == 1:
        return 'TNBC'

def add_hormone_info(ad, hormone_path):
    """
    Input: anndata file
    Output: anndata file with hormone information 
    """
    hormone = pd.read_csv(hormone_path, sep = ' ', header = 0)
    hormone['Classe'] = hormone.apply(hormone_classe, axis=1)
    sample_class_map = dict(zip(hormone['Sample_id'], hormone['Classe']))
    ad.obs['Hormone'] = [sample_class_map.get(sample_id, None) for sample_id in ad.obs.index]


#############################
######### Overlapp ##########
#############################


def adjust_bed_position(start, end, seq_len = 1344):
    """
        Adjust the position of the start and end in the BED file to make it 1344 bp long

        Parameters:
        - start: start of the DNA position
        - end: end of the DNA position

        Output:
        - Adjusted positions
    """

    mid = (start + end) // 2
    seq_start = mid - seq_len // 2
    seq_end = seq_start + seq_len
    if seq_start < 0:
        seq_start = 0
    return seq_start, seq_end


def new_bed(bed_file, output_file=f'./intermediate_file.bed'):
    """
        Create a new BED file with adjusted positions

        Parameters:
        - bed_file: BED file
        - output_file: intermediate BED file name

        Output:
        - BED file with adjusted positions
    """

    a = pybedtools.BedTool(bed_file)
    with open(output_file, 'w') as int_bed: 
        ### Write new intermediate bedfile
        for position in a:
            new_start, new_end = adjust_bed_position(int(position.start), int(position.end))
            int_bed.write(f"{position[0]}\t{new_start}\t{new_end}\n") 
    return output_file  # Returns the output file name


def global_overlap_distribution(data, percentage):
    """
        Create a histogram of the base pair overlapping in the BED file

        Parameters:
        - data: list of all overlapping DNA positions
        - percentage: percentage of overlapping

        Output:
        - PNG file of the histogram plot
    """

    plt.hist(data, bins=200, color='skyblue', edgecolor='black') 
    plt.xlabel('Number of Base Pair Overlap')
    plt.ylabel('Frequency')
    plt.title(f'Overlap Distribution ({percentage}%)')
    plt.xticks(range(0, 1345, 100)) 
    plt.savefig('./global_distribution.png', format='png')


def natural_sort_key(s):
    """
        Sort the dictionary keys by numerical order

        Parameters:
        - s: dictionary with keys to sort

        Output:
        - Ordered dictionary keys
    """

    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]


def chr_overlap_distribution(overlap_chr, percentage):
    """
        Create a density plot for the base pair overlap distribution from a BED file for each chromosome.

        Parameters:
        - overlap_chr: Dictionary with 'chrom' as keys and [list of overlaps] as values.
        - percentage: percentage of global overlap

        Output:
        - Density plot for each chromosome's base pair overlap.
    """

    plt.figure(figsize=(10, 6))
    palette = sns.color_palette("husl", len(overlap_chr))

    sorted_chromosomes = sorted(overlap_chr.keys(), key=natural_sort_key)
    
    # Loop through each chromosome and its corresponding values in the dictionary
    legend_elements = []  
    for i, chrom in enumerate(sorted_chromosomes):
        values = overlap_chr[chrom]
        # Convert values to numeric data
        values_numeric = pd.to_numeric(values)
        
        # Plot the density curve with a unique color for each chromosome
        sns.histplot(values_numeric, kde=True, color=palette[i], element='poly', 
                     fill=False, linewidth=2, alpha=0, label=chrom)
        # Create a legend element for each chromosome
        legend_elements.append(Line2D([0], [0], color=palette[i], lw=2, label=chrom))

    # Create the legend from the custom legend elements
    plt.legend(handles=legend_elements, title='Chromosome')

    plt.xlabel('Number of Base Pair Overlap')
    plt.ylabel('Frequency')
    plt.title(f'Overlap Distribution ({percentage}%) per Chromosome')
    plt.xticks(range(0, 1345, 100)) 
    plt.savefig('./chr_distribution.png', format='png')


def calculate_overlap_percentage(bed_file):
    """
        Calculate the percentage of overlap in the given BED file

        Parameters:
        - bed_file: BED file

        Output:
        - Percentage of overlap in the BED file
    """

    a = pybedtools.BedTool(bed_file)
    
    # Calculate overlap between features in the BED file
    overlaps = a.intersect(a, wo=True)
    
    # Set to store unique regions with overlap
    regions_with_overlap = set()

    # list to store the bp overlap
    overlap_bp = []
    overlap_chr = {}

    # Iterate over overlap results and add unique regions to the set
    for overlap in overlaps:
        chrom, start1, end1, _, start2, end2, overlap_size = overlap.fields[0:7]
        if start1 != start2 or end1 != end2:  # Make sure it's not the same region overlapping itself
            regions_with_overlap.add((chrom, start1, end1))
            regions_with_overlap.add((chrom, start2, end2))
            if chrom in overlap_chr:   
                overlap_chr[chrom].append(overlap_size)
                overlap_chr[chrom].sort(key=int)
            else:
                overlap_chr[chrom] = [overlap_size]
            overlap_bp.append(overlap_size)
    
    overlap_bp = sorted(overlap_bp,  key=int)

    # Calculate the percentage of regions that have overlap
    num_regions_with_overlap = len(regions_with_overlap)
    total_regions = len(a)
    overlap_percentage = (num_regions_with_overlap / total_regions) * 100 if total_regions > 0 else 0
    print(f'% of overlapping: {round(overlap_percentage, 2)}')
    return overlap_percentage, overlap_bp, overlap_chr
