
# Installation

To install this package run the following.
```
git clone https://github.com/AndrewLutsky/SCKDE.git
conda create -n sckde
pip install -e .
```

To use this package on any anndata object:
```
from sckde.sckde import sckde
import scanpy as sc

adata.obs['Density'] = sckde(adata, key='CD4')
sc.pl.umap(adata.obs, color = 'Density', title='CD4 Density')
```



# TODO
- Integration Tests
- Finish README
- Unit Tests
- Create plotting module
- Play around with bandwidths.
- Rewrite sckde arguments
    - allow for different bw/KDEpy methods
    - check default bw method for fit
    - allow for other KDEpy kernels
- Write docs
- Run experiments
    - How do other kernels do with weighted
    count data?
    - How does the package compare to nebulosa?
    - how does scipy compare to?
    - What does pseudotime look like along path?
