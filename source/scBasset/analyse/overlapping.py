import argparse
import warnings
from analyse_function import calculate_overlap_percentage, global_overlap_distribution, chr_overlap_distribution

warnings.filterwarnings('ignore')

"""
Script that calculates the percentage of overlap between chromosome positions (extended to 1344 bp) within a BED file. 
Possibility to visualize the global overlap distribution through a histogram as well as a chromosome-by-chromosome visualization. 
"""


def overlapp(bed_file, plot = False):
    overlap_percentage, overlap_bp, overlap_chr = calculate_overlap_percentage(bed_file)
    if plot == True:
        global_overlap_distribution(overlap_bp, round(overlap_percentage, 2))
        chr_overlap_distribution(overlap_chr,round(overlap_percentage, 2))

def main():
    parser = argparse.ArgumentParser(description='Calculate the overlap between DNA positions in a BED file')
    parser.add_argument('bed_path', type=str, help='Path to the BED file containing the DNA positions.')
    parser.add_argument('--plot', action='store_true', help='Specify if a plot is needed.')

    args = parser.parse_args()

    overlapp(bed_file=args.bed_path, plot=args.plot)
    
   
    
if __name__ == "__main__":
    main()
