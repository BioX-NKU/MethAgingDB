import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

base_dir = './'
output_dir = './output/'
os.makedirs(output_dir,exist_ok =True)

GSE_numbers = os.listdir(base_dir)
if '.ipynb_checkpoints' in GSE_numbers:
    GSE_numbers.remove('.ipynb_checkpoints')
error_GSE_numbers = []

for GSE_number in GSE_numbers:
        print(GSE_number)
        matrix = pd.read_csv(base_dir + f'/{GSE_number}/{GSE_number}_count_matrix.csv',index_col=0,sep='\t')
        matrix.fillna(0, inplace=True)
        # Perform PCA to reduce dimensions
        n_components = min(50,matrix.shape[1])
        pca = PCA(n_components=n_components)  # Reduce to 50 dimensions if the original dimensionality is very high
        pca_results = pca.fit_transform(matrix.T)  # Transpose df to have samples as rows

        pca_top2 = pca_results[:, :2]

        pca_df = pd.DataFrame(pca_top2, index=matrix.columns, columns=['PCA-1', 'PCA-2'])
        pca_df.to_csv(output_dir + f'/{GSE_number}_PCA_results.csv')
    