import pybedtools
import sys
import os 
import matplotlib.pyplot as plt


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


def display_distribution(data, percentage):
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

    # Iterate over overlap results and add unique regions to the set
    for overlap in overlaps:
        chrom, start1, end1, _, start2, end2, _ = overlap.fields[0:7]
        if start1 != start2 or end1 != end2:  # Make sure it's not the same region overlapping itself
            regions_with_overlap.add((chrom, start1, end1))
            regions_with_overlap.add((chrom, start2, end2))
            overlap_bp.append(overlap.fields[6])
    
    overlap_bp = sorted(overlap_bp,  key=int)

    # Calculate the percentage of regions that have overlap
    num_regions_with_overlap = len(regions_with_overlap)
    total_regions = len(a)

    overlap_percentage = (num_regions_with_overlap / total_regions) * 100 if total_regions > 0 else 0

    display_distribution(overlap_bp, round(overlap_percentage, 2))
    
    print(f'% of overlapping: {round(overlap_percentage, 2)}')
    
    
    
if __name__ == "__main__":
    bed_file = sys.argv[1]
    new_bed_file = new_bed(bed_file)
    calculate_overlap_percentage(new_bed_file)
    os.remove(new_bed_file)
