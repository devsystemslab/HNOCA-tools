from typing import Optional
import anndata as ad
import pandas as pd

from .auroc import annotate_snap

from .utils import get_markers, get_annot_df


def annotate_hierarchy(
    adata: ad.AnnData,
    marker_hierarchy: dict,
    group_name: str,
    layer: Optional[str] = None,
    min_expr: float = 0.1,
    **kwargs
) -> dict:
    """
    Annotate clusters based on a manually defined cell type and marker hierarchy.

    Args:
        adata: AnnData object
        marker_hierarchy: dict arker genes for each celltype arranged
            hierarchically.
        group_name: Name of the column in adata.obs that contains the cluster labels
        layer: Layer in adata to use for expression
        **kwargs: Additional arguments to pass to the annotation function.

    Returns:
        Dict with assignments and metrics
    """

    # Annotate at each level of the hierarchy
    assignment_hierarchy = annotate_levels(
        adata, marker_hierarchy, group_name, **kwargs
    )

    return dict(
        assignments=get_annot_df(assignment_hierarchy, group_name, min_expr=min_expr),
        metrics=assignment_hierarchy,
    )


def annotate_levels(
    adata,
    marker_hierarchy,
    group_name,
    level=0,
    assignment_levels=None,
    layer=None,
    **kwargs
):
    """Recursively annotatates all levels of a marker hierarchy."""
    level += 1
    level_name = "level_" + str(level)
    marker_dict = get_markers(marker_hierarchy)
    assignments = annotate(adata, marker_dict, group_name, layer=layer, **kwargs)

    if assignment_levels is None:
        assignment_levels = {}

    if level_name not in assignment_levels.keys():
        assignment_levels[level_name] = pd.DataFrame()

    assignment_levels[level_name] = pd.concat(
        [assignment_levels[level_name], assignments], axis=0
    )

    for subtype in assignments["class"].unique():

        if "subtypes" not in marker_hierarchy[subtype].keys():
            continue

        # Subset adata
        subtype_groups = assignments[group_name][
            assignments["class"] == subtype
        ].astype(str)
        subtype_adata = adata[adata.obs[group_name].isin(subtype_groups)]

        # Recursively annotate
        assignment_levels = annotate_levels(
            subtype_adata,
            marker_hierarchy[subtype]["subtypes"],
            group_name,
            level=level,
            assignment_levels=assignment_levels,
            layer=layer,
        )

    return assignment_levels


def annotate(
    adata: ad.AnnData,
    marker_dict: dict,
    group_name: str,
    layer: Optional[str] = None,
    **kwargs
) -> pd.DataFrame:
    """
    Annotate clusters based on a manually defined cell type markers.

    Args:
        adata: AnnData object
        marker_dict: Dict with marker genes for each celltype
        group_name: Name of the column in adata.obs that contains the cluster labels
        layer: Layer in adata to use for expression
        **kwargs: Additional arguments to pass to the annotation function.

    Returns:
        pd.DataFrame with assignments
    """
    assignments = annotate_snap(adata, marker_dict, group_name, layer=layer, **kwargs)
    # Join cluster-level results with adata
    assignments = assignments.reset_index(names=group_name)
    return assignments
