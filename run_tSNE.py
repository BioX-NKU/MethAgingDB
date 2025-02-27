base_dir = './'
gse_id = 'GSE14788'
matrix = pd.read_csv(base_dir + f'/{gse_id}_count.csv',index_col=0,sep='\t')
meta = pd.read_csv(base_dir + f'/{gse_id}_meta.csv',index_col=0,sep='\t')
matrix.fillna(0, inplace=True)
# Perform PCA to reduce dimensions
n_components = min(50,matrix.shape[1])
pca = PCA(n_components=n_components)  # Reduce to 50 dimensions if the original dimensionality is very high
pca_results = pca.fit_transform(matrix.T)  # Transpose df to have samples as rows

perplexity_value = min(30, matrix.shape[1] - 1) 

tsne = TSNE(n_components=2, random_state=42)
tsne_results = tsne.fit_transform(pca_results) 
tsne_df = pd.DataFrame(tsne_results, index=matrix.columns, columns=['tSNE-1', 'tSNE-2'])
tsne_df = tsne_df.join(meta)
tsne_df.to_csv(base_dir +f'/{gse_id}_tSNE_results.csv')

