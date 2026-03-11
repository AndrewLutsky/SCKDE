import scanpy as sc
from sckde.sckde import sckde

pbmc = sc.read_h5ad("/Users/andrewlutsky/Desktop/W\
ork/sckde/tests/test_data/pbmc.h5ad")
pbmc.obs['density'] = sckde(pbmc, keys=["CD8A", "CCR7"])
sc.pl.umap(pbmc, color='density')

