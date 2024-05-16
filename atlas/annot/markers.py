import pandas as pd

import jax
from jax import numpy as jnp

from functools import partial

from .utils import get_expr, matrix_to_long_df
from .auroc import auc_expr


def find_markers(adata, group_name, features=None, layer=None):
    """
    Find markers for each cluster.

    Parameters
    ----------
    adata
        AnnData object
    group_name
        Name of the column in adata.obs that contains the cluster labels
    features
        List of features to use for marker identification
    layer
        Layer in adata to use for expression

    Returns
    -------
    DataFrame with AUROC and detection ratio for each gene and cluster.
    """
    # Get expression matrix
    expr, features = get_expr(adata, features=features, layer=layer)
    # Compute AUROC and fraction nonzero for marker features
    metrics = auc_expr(
        adata,
        group_name,
        features=features,
        compute_frac_nonzero_out=True,
        apply_fun="numpy",
    )
    # Make dataframes and join
    expr_df = matrix_to_long_df(
        metrics["auroc"], features=metrics["features"], groups=metrics["groups"]
    ).rename(columns={"value": "auroc"})

    auc_df = matrix_to_long_df(
        metrics["frac_nonzero"], features=metrics["features"], groups=metrics["groups"]
    ).rename(columns={"value": "frac_nonzero"})

    out_expr_df = matrix_to_long_df(
        metrics["frac_nonzero_out"],
        features=metrics["features"],
        groups=metrics["groups"],
    ).rename(columns={"value": "frac_nonzero_out"})

    marker_df = expr_df.merge(auc_df, on=["group", "feature"]).merge(
        out_expr_df, on=["group", "feature"]
    )
    # Compute detection ratio
    marker_df["detection_ratio"] = (
        marker_df["frac_nonzero"] / marker_df["frac_nonzero_out"]
    )

    return marker_df.sort_values(["group", "auroc"], ascending=False)
