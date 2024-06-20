import matplotlib.pyplot as plt
import anndata
import pandas as pd

def plot_probability_distribution(file_path, output_path, bins=50, title='Probability Distribution', xlabel='Probabilities', ylabel='Density'):
    # Load probabilities from the text file while skipping the first line
    probabilities = []
    with open(file_path, 'r') as f:
        next(f)  # Skip the first line
        for line in f:
            probabilities.extend(map(float, line.strip().split()))

    # Create the histogram
    plt.hist(probabilities, bins=bins, density=True, alpha=0.7, color='blue')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.savefig(output_path)  # Save the plot to a file
    plt.close()  # Close the plot to release resources

# Spécifiez le chemin de sortie souhaité pour enregistrer le tracé
output_file_path = 'probability_distribution_tutorial_train_epoch1000_plot.png'

# Appel de la fonction pour générer et enregistrer le tracé
plot_probability_distribution('probabilities_tutorial_train_epoch1000_half.txt', output_file_path)

