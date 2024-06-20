import pandas as pd
import umap
import matplotlib.pyplot as plt
import sys
import numpy as np

def visualize_umap(csv_file):
    # Load CSV data using pandas
    data = pd.read_csv(csv_file, index_col=0)  # Use the first column as index

    # Transpose the data so that columns become rows and vice versa
    data_transposed = data.transpose()

    # Split the data into two sets based on columns
    cols_last_10 = data_transposed.iloc[:, -10:].values
    cols_rest = data_transposed.iloc[:, :-10].values

    # Dimension reduction with UMAP
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(data_transposed)

    # Plot UMAP
    plt.scatter(embedding[:, 0], embedding[:, 1], c=np.array([0]*len(cols_rest) + [1]*len(cols_last_10)), cmap=plt.cm.Set1)
    plt.colorbar(boundaries=[0, 1], ticks=[0, 1], label='Color')
    plt.title('UMAP Visualization')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.show()

# Example of usage
if __name__ == "__main__":
    csv_file = sys.argv[1]
    visualize_umap(csv_file)

