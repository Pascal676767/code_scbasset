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


def embedding(model, results_path):
    """
        Extract the weights of the final dense layer

        Parameters:
        - model: trained model
        - results_path: path to save the results

        Output:
        - CSV projection file
        
    """

    proj = get_cell_embedding(model) # get_cell_embedding function
    pd.DataFrame(proj).to_csv(f'{results_path}/projection_cage.csv')


def plotting_embedding(ad_file, color_by, results_path):
    """
        Visualize the embedding of the last dense layer weights, colored by a specified attribute

        Parameters:
        - ad_file: file with samples and features
        - color_by: column name to color by
        - results_path: path where the PNG will be saved

        Output:
        - PNG plot of the embeddings 
    """

    # Load the projection
    projection_path = f'{results_path}/projection_cage.csv'
    ad_file.obsm['projection'] = pd.read_csv(projection_path, index_col=0).values
    
    # Compute neighbors and UMAP
    sc.pp.neighbors(ad_file, use_rep='projection')
    sc.tl.umap(ad_file)
    sc.tl.leiden(ad_file)
    
    # Create a figure and axis
    f, ax = plt.subplots(figsize=(10, 10))
    
    # Plot UMAP embedding, disable immediate display
    sc.pl.umap(ad_file, color=color_by, ax=ax, size=350, title='embedding', show=False)
    
    # Adjust the layout
    plt.tight_layout()

    output_file = f'{results_path}/embedding_{color_by}.png'
    
    # Save the plot and handle potential exceptions
    try:
        plt.savefig(output_file, format='png')
        print(f"Embedding plot saved to {output_file}")
    except Exception as e:
        print(f"Failed to save embedding plot: {e}")

    plt.close()
    

def add_patient_info(ad):
    """
        Add patient or cell line information to every sample in the AnnData file

        Parameters:
        - ad: AnnData file with samples and features

        Output:
        - Updated AnnData file with patient and cell line information
    """

    ad.obs['Patient'] = np.where(ad.obs.index.str.startswith('s'), 'patient', 'cell_line')
    ad.obs['Patient'] = ad.obs['Patient'].astype('category')

def hormone_classe(row):
    """
        Add hormone status information to every sample corresponding to a patient in the AnnData file

        Parameters:
        - row: single row of the AnnData file containing sample information

        Output:
        - Subgroup name of the breast cancer type
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
        Add hormone status information to every sample corresponding to a patient in the AnnData file

        Parameters:
        - ad: AnnData file with samples and features 

        Output:
        - Updated AnnData file with hormonal status information for every patient
    """

    hormone = pd.read_csv(hormone_path, sep = '\t', header = 0)
    hormone['Classe'] = hormone.apply(hormone_classe, axis=1)
    sample_class_map = dict(zip(hormone['Sample_id'], hormone['Classe']))
    ad.obs['Hormone'] = [sample_class_map.get(sample_id, None) for sample_id in ad.obs.index]


################################
######### Overlapping ##########
################################


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
        - global overlapping
        - overlapping per chromosomes
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
