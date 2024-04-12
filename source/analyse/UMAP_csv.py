import pandas as pd
import umap
import matplotlib.pyplot as plt
import sys
import numpy as np

def visualize_umap(csv_file):
    # Charger les données CSV en utilisant pandas
    data = pd.read_csv(csv_file, index_col=0)  # Utilisez la première colonne comme index

    # Transposer les données pour que les colonnes deviennent des lignes et vice versa
    data_transposed = data.transpose()

    # Séparer les données en deux ensembles en fonction des colonnes
    cols_last_10 = data_transposed.iloc[:, -10:].values
    cols_rest = data_transposed.iloc[:, :-10].values

    # Réduction de dimension avec UMAP
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(data_transposed)

    # Plot de l'UMAP
    plt.scatter(embedding[:, 0], embedding[:, 1], c=np.array([0]*len(cols_rest) + [1]*len(cols_last_10)), cmap=plt.cm.Set1)
    plt.colorbar(boundaries=[0, 1], ticks=[0, 1], label='Color')
    plt.title('UMAP Visualization')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.show()
# Exemple d'utilisation
if __name__ == "__main__":
    csv_file = sys.argv[1]
    visualize_umap(csv_file)
