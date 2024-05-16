# Human Neural Organoid Cell Atlas Toolbox
#### 🛠️ The Swiss Army Knive of the Single Cell Cartographer

This package provides a set of tools we used to generate and analyze the Human Neural Organoid Cell Atlas. Among other things, it provides functions to:

- Rapidly annotate cell types based on marker genes
- Map query data to the reference atlas
- Transfer annotations between datasets
- Compute 'presence scores' for query data based on the reference atlas
- Perform differential expression analysis


## Quick start

### 🖋️ Annotation 

We developed [snapseed](https://github.com/devsystemslab/snapseed) to rapidly annotate the HNOCA. It annotates cells based on manually defined sets of marker genes for individual cell types or cell type hierarchies. It is fast (i.e. GPU-accelerated) and simple to enable annotation of very large datasets.

```python
import hnoca.snapseed as snap

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

### 🗺️ Mapping

For reference mapping, we mostly rely on [scPoli](https://docs.scarches.org/en/latest/scpoli_surgery_pipeline.html) and [scANVI](https://docs.scvi-tools.org/en/1.1.1/user_guide/models/scanvi.html). Based on pretrained models, we here provide a simple interface to map query data to the reference atlas.

```python