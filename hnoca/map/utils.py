import warnings
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad


def prepare_features(query_adata, ref_model):
    """
    Prepare the features of the query dataset to match the reference dataset
    """
    ref_features = ref_model.adata.var_names
    query_features = query_adata.var_names
    missing_features = np.setdiff1d(ref_features, query_features)

    if len(missing_features) == 0:
        return query_adata

    print(
        f"Warning: Query dataset is missing {len(missing_features)} features from the reference dataset. Adding missing features as zero-filled columns."
    )

    new_adata = query_adata.copy()
    new_adata = new_adata[:, np.isin(query_features, ref_features)]

    extra_X = sparse.csr_matrix(
        (new_adata.shape[0], len(missing_features)), dtype=new_adata.X.dtype
    )
    new_X = sparse.hstack((new_adata.X, extra_X))

    new_layers = {}
    for layer in new_adata.layers.keys():
        new_layers[layer] = sparse.hstack((new_adata.layers[layer], extra_X))

    new_var = pd.DataFrame(
        {"symbol": list(new_adata.var_names) + list(missing_features)}
    ).set_axis(list(new_adata.var_names) + list(missing_features))

    new_adata = ad.AnnData(
        X=new_X,
        obs=new_adata.obs.copy(),
        var=new_var,
        obsm=new_adata.obsm.copy(),
        layers=new_layers,
    )
    new_adata = new_adata[:, ref_features].copy()

    return new_adata
