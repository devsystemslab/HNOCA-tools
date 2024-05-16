import pandas as pd

from .trinarize import annotate_cytograph
from .auroc import annotate_snap

from .utils import get_markers, get_annot_df


def annotate_hierarchy(
    adata,
    marker_hierarchy,
    group_name,
    method="auroc",
    layer=None,
    min_expr=0.1,
    **kwargs
):
    """
    Annotate clusters based on a manually defined cell type and marker hierarchy.

    Parameters
    ----------
    adata
        AnnData object
    marker_hierarchy
        Dict with marker genes for each celltype arranged hierarchically.
    group_name
        Name of the column in adata.obs that contains the cluster labels
    method
        Method to use for annotation. Options are "auroc" and "trinatize".
    layer
        Layer in adata to use for expression
    **kwargs
        Additional arguments to pass to the annotation function.
    """

    # Annotate at each level of the hierarchy
    assignment_hierarchy = annotate_levels(
        adata, marker_hierarchy, group_name, method=method, **kwargs
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
    method="auroc",
    layer=None,
    **kwargs
):
    """Recursively annotatates all levels of a marker hierarchy."""
    level += 1
    level_name = "level_" + str(level)
    marker_dict = get_markers(marker_hierarchy)
    assignments = annotate(
        adata, marker_dict, group_name, method=method, layer=layer, **kwargs
    )

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
            method=method,
            layer=layer,
        )

    return assignment_levels


def annotate(adata, marker_dict, group_name, layer=None, **kwargs):
    """
    Annotate clusters based on a manually defined cell type markers.

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
    **kwargs
        Additional arguments to pass to the annotation function.
    """
    assignments = annotate_snap(adata, marker_dict, group_name, layer=layer, **kwargs)
    # Join cluster-level results with adata
    assignments = assignments.reset_index(names=group_name)
    return assignments
