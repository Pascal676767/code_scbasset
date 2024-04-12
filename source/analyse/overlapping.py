import pybedtools
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D
import re


def adjust_bed_position(start, end, seq_len = 1344):
    """
    Adjust the position of the start and end in the BED file to make it 1344 bp long
    Input: A position with a start and an end
    Output: Adjusted position
    """
    mid = (start + end) // 2
    seq_start = mid - seq_len // 2
    seq_end = seq_start + seq_len
    if seq_start < 0:
        seq_start = 0
    return seq_start, seq_end

def new_bed(bed_file, output_file=f'int_{sys.argv[1]}'):
    """
    Create a new BED file with adjusted positions
    Input: BED file
    Output: New BED file
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
    Input: List of sorted overlap values
    Output: Histogram
    """
    plt.hist(data, bins=200, color='skyblue', edgecolor='black') # bins determines the number of bars in the histogram
    plt.xlabel('Number of Base Pair Overlap')
    plt.ylabel('Frequency')
    plt.title(f'Overlap Distribution ({percentage}%)')
    plt.xticks(range(0, 1345, 100)) 
    plt.show()

def natural_sort_key(s):
    """Clé de tri naturel pour trier les chromosomes par ordre numérique."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]


def chr_overlap_distribution(overlap_chr, percentage):
    """
    Create a density plot for the base pair overlap distribution from a BED file for each chromosome.
    Input: Dictionary with 'chrom' as keys and [list of values] as values.
    Output: Density plot for each chromosome's base pair overlap.
    """
    plt.figure(figsize=(10, 6))
    palette = sns.color_palette("husl", len(overlap_chr))

    sorted_chromosomes = sorted(overlap_chr.keys(), key=natural_sort_key)
    
    # Loop through each chromosome and its corresponding values in the dictionary
    legend_elements = []  # list to hold the legend elements
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
    plt.show()


def calculate_overlap_percentage(bed_file):
    """
    Calculate the percentage of overlap in the given BED file
    Input: BED file
    Output: Percentage
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
    print(overlap_bp)
    overlap_percentage = (num_regions_with_overlap / total_regions) * 100 if total_regions > 0 else 0

    global_overlap_distribution(overlap_bp, round(overlap_percentage, 2))
    chr_overlap_distribution(overlap_chr,round(overlap_percentage, 2) )
    
    print(f'% of overlapping: {round(overlap_percentage, 2)}')
    
    
    
if __name__ == "__main__":
    bed_file = sys.argv[1]
    new_bed_file = new_bed(bed_file)
    calculate_overlap_percentage(new_bed_file)
    os.remove(new_bed_file)
