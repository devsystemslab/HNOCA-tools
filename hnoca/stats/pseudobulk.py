from typing import Optional
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from decoupler import get_pseudobulk


def create_pseudobulk(
    adata: ad.AnnData,
    sample_key: str = "batch",
    group_key: Optional[str] = None,
    layer: Optional[str] = None,
    mode: str = "sum",
    min_cells: int = 10,
    min_counts: int = 1000,
    **kwargs,
) -> ad.AnnData:
    """
    Create pseudobulk data from anndata object

    Args:
        adata: AnnData object
        sample_key: Column name in adata.obs that contains the sample ID
        group_key: Column name in adata.obs that contains the group ID
        layer: Layer to use for pseudobulk data. If None, use `adata.X`
        mode: Method to aggregate data. Default is 'sum'.
        min_cells: Filter to remove samples by a minimum number of cells in a sample-group pair.
        min_counts: Filter to remove samples by a minimum number of summed counts in a sample-group pair.
        **kwargs: Additional arguments to pass to `decoupler.get_pseudobulk()`


    Returns:
        AnnData object with pseudobulk data
    """
    if layer is not None:
        adata.X = adata.layers[layer].copy()

    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A.ravel()

    if group_key is not None:
        adata.obs["group"] = adata.obs[group_key]
    else:
        adata.obs["group"] = "group"

    adata_pb = get_pseudobulk(
        adata,
        sample_col=sample_key,
        groups_col="group",
        layer=layer,
        mode=mode,
        min_cells=min_cells,
        min_counts=min_counts,
        **kwargs,
    )
    adata_pb.X = np.round(adata_pb.X)
    adata_pb.obs["n_genes_median"] = [
        adata.obs.groupby("batch").median()["n_genes"].loc[i.split("_onegroup")[0]]
        for i in adata_pb.obs.index
    ]
    adata_pb.obs["n_genes_std"] = [
        adata.obs.groupby("batch").std()["n_genes"].loc[i.split("_onegroup")[0]]
        for i in adata_pb.obs.index
    ]

    return adata_pb
