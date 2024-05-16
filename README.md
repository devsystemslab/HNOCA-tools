# Human Neural Organoid Cell Atlas Toolbox
#### The Swiss Army Knive of the Single Cell Cartographer

This package provides a set of tools we used to generate and analyze the Human Neural Organoid Cell Atlas. Among other things, it provides functions to:

- Annotate cell types based on marker genes
- Map query data to the reference atlas
- Transfer annotations between datasets
- Compute 'presence scores' for query data based on the reference atlas
- Perform differential expression analysis


## Quick start

### Annotation



```python
import snapseed as snap
from snapseed.utils import read_yaml

# Read in the marker genes
marker_genes = read_yaml("marker_genes.yaml")

# Annotate anndata objects
snap.annotate(
    adata,
    marker_genes,
    group_name="clusters",
    layer="lognorm",
)

# Or for more complex hierarchies
snap.annotate_hierarchy(
    adata,
    marker_genes,
    group_name="clusters",
    layer="lognorm",
)
```
