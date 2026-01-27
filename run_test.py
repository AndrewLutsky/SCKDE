import scanpy as sc
from sckde import _sckde_multi
import seaborn as sb 
import matplotlib.pyplot as plt 
adata = sc.read_h5ad("tests/test.h5ad")
adata.var_names_make_unique()
print(adata.var_names[0:100])
# Filled contour
adata.obs['Density'] = _sckde_multi(adata, ['CDK16'])
#cnt = ax.contourf(df['x'], df['y'], df['Density'])
print(adata.obs)
sc.pl.umap(adata, color='Density', size=100, alpha=0.5)
sc.pl.umap(adata, color='CDK16', size=100, alpha=0.5)
plt.show()
