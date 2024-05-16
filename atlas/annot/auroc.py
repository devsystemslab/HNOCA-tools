import pandas as pd
import numpy as np

import jax
from jax import numpy as jnp

from functools import partial
from sklearn import preprocessing

from .utils import (
    dict_to_binary,
    get_expr,
    frac_nonzero,
    masked_max,
    masked_mean,
    jit_frac_nonzero,
)


def annotate_snap(
    adata,
    marker_dict,
    group_name,
    layer=None,
    auc_weight=0.5,
    expr_weight=0.5,
    marker_summary_fun="max",
    apply_fun="vmap",
):
    """
    Annotate cell types based on AUROC and expression of predefined marker genes.

    Parameters
    ----------
    adata
        AnnData object
    marker_dict
        Dict with marker genes for each celltype
    group_name
        Name of the column in adata.obs that contains the cluster labels
    layer
        Layer in adata to use for expression
    auc_weight
        Weight to give to AUROC in the final score
    expr_weight
        Weight to give to detection rate in the final score
    marker_summary_fun
        Function to use to summarize over markers for the same annotation. Options are "max" and "mean".
    apply_fun
        Function to use for applying the auroc function. Options are "vmap" and "numpy". vmap is a lot faster but can cause OOM errors for big matrices.
    """
    # Reformat marker_dict into binary matrix
    marker_mat = dict_to_binary(marker_dict)

    # Compute AUROC and fraction nonzero for marker features
    features = marker_mat.columns
    metrics = auc_expr(
        adata, group_name, features=features, layer=layer, apply_fun=apply_fun
    )

    marker_mat = marker_mat.loc[:, metrics["features"]]
    if marker_summary_fun == "max":
        auc_max = masked_max(metrics["auroc"], marker_mat.values)
        expr_max = masked_max(metrics["frac_nonzero"], marker_mat.values)
    elif marker_summary_fun == "mean":
        auc_max = masked_mean(metrics["auroc"], marker_mat.values)
        expr_max = masked_mean(metrics["frac_nonzero"], marker_mat.values)
    else:
        raise ValueError(
            f"Invalid marker_summary_fun: {marker_summary_fun}. "
            "Current options are 'max' and 'mean'."
        )
    # Combine metrics
    assignment_scores = (auc_weight * auc_max + expr_weight * expr_max) / (
        auc_weight + expr_weight
    )
    assign_idx = jnp.argmax(assignment_scores, axis=0)
    # Mask out genes that are not expressed in any cell
    assign_class = marker_mat.index[assign_idx]

    assign_df = pd.DataFrame(
        {
            "class": assign_class,
            "score": assignment_scores[assign_idx, jnp.arange(auc_max.shape[1])],
            "auc": auc_max[assign_idx, jnp.arange(auc_max.shape[1])],
            "expr": expr_max[assign_idx, jnp.arange(expr_max.shape[1])],
        },
        index=metrics["groups"],
    )

    return assign_df


def auc_expr(
    adata,
    group_name,
    features=None,
    layer=None,
    compute_auroc=True,
    compute_frac_nonzero=True,
    compute_frac_nonzero_out=False,
    apply_fun="vmap",
):
    """Computes AUROC and fraction nonzero for each gene in an adata object."""
    # Turn string groups into integers
    le = preprocessing.LabelEncoder()
    le.fit(adata.obs[group_name])

    # Compute AUROC and fraction nonzero
    groups = jnp.array(le.transform(adata.obs[group_name]))
    expr, features = get_expr(adata, features=features, layer=layer)
    auroc, frac_nonzero, frac_nonzero_out = expr_auroc_over_groups(
        expr,
        groups,
        compute_auroc=compute_auroc,
        compute_frac_nz=compute_frac_nonzero,
        compute_frac_nz_out=compute_frac_nonzero_out,
        apply_fun=apply_fun,
    )

    return dict(
        frac_nonzero=frac_nonzero,
        frac_nonzero_out=frac_nonzero_out,
        auroc=auroc,
        features=features,
        groups=le.classes_,
    )


@jax.jit
def jit_auroc(x, groups):
    # TODO: compute frac nonzero here to avoid iterating twice

    # sort scores and corresponding truth values
    desc_score_indices = jnp.argsort(x)[::-1]
    x = x[desc_score_indices]
    groups = groups[desc_score_indices]

    # x typically has many tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate a value for the end of the curve.
    distinct_value_indices = jnp.array(jnp.diff(x) != 0, dtype=jnp.int32)
    threshold_mask = jnp.r_[distinct_value_indices, 1]

    # accumulate the true positives with decreasing threshold
    tps_ = jnp.cumsum(groups)
    fps_ = 1 + jnp.arange(groups.size) - tps_

    # mask out the values that are not distinct
    tps = jnp.sort(tps_ * threshold_mask)
    fps = jnp.sort(fps_ * threshold_mask)
    tps = jnp.r_[0, tps]
    fps = jnp.r_[0, fps]
    fpr = fps / fps[-1]
    tpr = tps / tps[-1]
    area = jnp.trapz(tpr, fpr)
    return area


vmap_auroc = jax.vmap(jit_auroc, in_axes=[1, None])


def numpy_auroc(x, groups):
    jit_auroc_ = partial(jit_auroc, groups=groups)
    auroc = np.apply_along_axis(jit_auroc_, axis=0, arr=x)
    return auroc


def expr_auroc_over_groups(
    expr,
    groups,
    compute_auroc=True,
    compute_frac_nz=True,
    compute_frac_nz_out=False,
    apply_fun="vmap",
):
    """Computes AUROC for each group separately."""
    auroc = np.zeros((groups.max() + 1, expr.shape[1]))
    frac_nz = np.zeros((groups.max() + 1, expr.shape[1]))
    frac_nz_out = np.zeros((groups.max() + 1, expr.shape[1]))

    if apply_fun == "vmap":
        auroc_fun = vmap_auroc
        fracnz_fun = jit_frac_nonzero
    elif apply_fun == "numpy":
        auroc_fun = numpy_auroc
        fracnz_fun = frac_nonzero

    for group in range(groups.max() + 1):
        if compute_auroc:
            auroc[group, :] = np.array(auroc_fun(expr, groups == group))
        if compute_frac_nz:
            frac_nz[group, :] = np.array(fracnz_fun(expr[groups == group, :]))
        if compute_frac_nz_out:
            frac_nz_out[group, :] = np.array(fracnz_fun(expr[groups != group, :]))

    return auroc, frac_nz, frac_nz_out
