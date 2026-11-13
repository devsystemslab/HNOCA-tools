import anndata as ad
import numpy as np
import pandas as pd

from hnoca.utils import check_deps


def create_pseudobulk(
    adata: ad.AnnData,
    sample_key: str = "batch",
    group_key: str | None = None,
    layer: str | None = None,
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
        **kwargs: Additional arguments to pass to `decoupler.pp.pseudobulk()`


    Returns
    -------
        AnnData object with pseudobulk data
    """
    check_deps("decoupler")
    import decoupler as dc

    adata_pb = dc.pp.pseudobulk(adata, sample_col=sample_key, groups_col=group_key, layer=layer, mode=mode, **kwargs)
    # every pseudobulk sample has at least one cell; decoupler 2 rejects min_cells < 1
    dc.pp.filter_samples(adata_pb, min_cells=max(min_cells, 1), min_counts=min_counts)
    adata_pb.X = np.round(adata_pb.X)

    X = adata.layers[layer] if layer is not None else adata.X
    n_genes = pd.Series(np.asarray((X > 0).sum(axis=1)).ravel(), index=adata.obs_names)
    n_genes_by_sample = n_genes.groupby(adata.obs[sample_key], observed=True)
    adata_pb.obs["n_genes_median"] = adata_pb.obs[sample_key].map(n_genes_by_sample.median()).to_numpy()
    adata_pb.obs["n_genes_std"] = adata_pb.obs[sample_key].map(n_genes_by_sample.std()).to_numpy()

    return adata_pb
