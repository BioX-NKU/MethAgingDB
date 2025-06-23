import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def filter_sex_chromosomes(df, cg_chromosome_mapping):
   
    # Ensure cgIDs are consistent between the data and mapping
    common_cgids = df.index.intersection(cg_chromosome_mapping.index)
    df = df.loc[common_cgids]
    cg_chromosome_mapping = cg_chromosome_mapping.loc[common_cgids]

    # Filter out sex chromosomes
    filtered_df = df.loc[~cg_chromosome_mapping.CHR.isin(['X', 'Y'])]
    return filtered_df
   

base_dir = './'
output_dir = './output/'
os.makedirs(output_dir,exist_ok =True)

dataset_all = pd.read_csv('./aging_datasets_basic_information.csv',encoding='latin-1',index_col=0)
GPL_df = pd.read_csv('./GPL8490.txt',sep='\t',skiprows=38)
GPL_df = GPL_df.rename(columns=lambda x: x.upper())
cg_chromosome_mapping = GPL_df[['ID','CHR']]
cg_chromosome_mapping = cg_chromosome_mapping.set_index('ID')


dataset_all = dataset_all[dataset_all.GPL == 'GPL8490']
GSE_numbers = dataset_all.GSE_number

if '.ipynb_checkpoints' in GSE_numbers:
    GSE_numbers.remove('.ipynb_checkpoints')
error_GSE_numbers = []

for GSE_number in GSE_numbers:
        print(GSE_number)
        matrix = pd.read_csv(base_dir + f'/{GSE_number}/{GSE_number}_count_matrix.csv',index_col=0,sep='\t')
        matrix.fillna(0, inplace=True)
        matrix = filter_sex_chromosomes(matrix, cg_chromosome_mapping)

        # Perform PCA to reduce dimensions
        n_components = min(50,matrix.shape[1])
        pca = PCA(n_components=n_components)  # Reduce to 50 dimensions if the original dimensionality is very high
        pca_results = pca.fit_transform(matrix.T)  # Transpose df to have samples as rows
     
        perplexity_value = min(30, matrix.shape[1] - 1)  # 选择一个默认值，但确保它小于样本数

        tsne = TSNE(n_components=2, random_state=42)
        tsne_results = tsne.fit_transform(pca_results) 
        tsne_df = pd.DataFrame(tsne_results, index=matrix.columns, columns=['tSNE-1', 'tSNE-2'])
        tsne_df.to_csv(output_dir + f'/{GSE_number}_tSNE_results.csv')
    