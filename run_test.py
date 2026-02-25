import scanpy as sc
from sckde import sckde 
import seaborn as sb 
import matplotlib.pyplot as plt 
adata = sc.read_h5ad("tests/test.h5ad")
adata.var_names_make_unique()
# Filled contour
print(sckde(adata, ['CDK16']))
adata.obs['Density'] = sckde(adata, ['CDK16'])
print(adata.obs)
sc.pl.umap(adata, color='Density', size=100, alpha=0.5)
#sc.tl.embedding_density(adata, basis='umap', key='CDK16')
#sc.pl.embedding_density(adata)
plt.show()
