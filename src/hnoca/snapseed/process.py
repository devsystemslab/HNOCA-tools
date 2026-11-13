from typing import Any

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.colors import hsv_to_rgb, rgb_to_hsv, to_hex, to_rgb

from hnoca._logging import logger


def _collect_marker_genes(node: dict[str, Any], gene_set: set[str] | None = None) -> set[str]:
    if gene_set is None:
        gene_set = set()
    if isinstance(node, dict):
        if "marker_genes" in node:
            gene_set.update(node["marker_genes"])
        if "subtypes" in node:
            for subnode in node["subtypes"].values():
                _collect_marker_genes(subnode, gene_set)
        # Also process any other nested dictionaries
        for key, value in node.items():
            if isinstance(value, dict) and key not in ["marker_genes", "subtypes"]:
                _collect_marker_genes(value, gene_set)
    return gene_set


def display_missing_genes(data: dict[str, Any], adata: ad.AnnData) -> None:
    """
    Log the marker genes of a hierarchical marker dict that are missing from `adata.var_names`.

    Args:
        data: Marker genes for each cell type, arranged hierarchically.
        adata: AnnData object.
    """
    all_genes = _collect_marker_genes(data)
    missing_genes = all_genes - set(adata.var_names)

    if missing_genes:
        logger.info("Genes not found in adata.var:")
        for gene in sorted(missing_genes):
            logger.info(gene)
    else:
        logger.info("All marker genes are present in adata.var.")


def derive_child_colors(
    parent_color: str,
    num_children: int,
    rotation_step: float,
    saturation_adjust: float,
    brightness_adjust: float,
    min_hue_separation: float,
    used_colors: set | None = None,
) -> list[str]:
    """
    Generate child colors by rotating the hue of the parent color and adjusting saturation and brightness.

    Args:
        parent_color: Parent color in hex format.
        num_children: Number of child colors to generate.
        rotation_step: Step size for hue rotation.
        saturation_adjust: Adjustment factor for saturation (positive to increase, negative to decrease).
        brightness_adjust: Adjustment factor for brightness (positive to increase, negative to decrease).
        min_hue_separation: Minimum separation between hue values of child colors and between parent and child hues.
        used_colors: Colors already used in the hierarchy, to avoid duplicates. Updated in place.

    Returns
    -------
        List of child colors in hex format.
    """
    if parent_color is None:
        raise ValueError("Parent color cannot be None when generating child colors.")
    if num_children <= 0:
        return []

    rotation_step = max(rotation_step, min_hue_separation)

    parent_hsv = rgb_to_hsv(np.array(to_rgb(parent_color)))
    parent_hue = parent_hsv[0]
    child_colors = []
    used_colors = used_colors if used_colors is not None else set()

    for i in range(num_children):
        hue = (parent_hue + i * rotation_step) % 1.0
        saturation = min(max(parent_hsv[1] + saturation_adjust * ((-1) ** i), 0.0), 1.0)
        brightness = min(max(parent_hsv[2] + brightness_adjust * ((-1) ** i), 0.0), 1.0)

        # Ensure minimum hue separation between parent and child
        if abs(hue - parent_hue) < min_hue_separation:
            hue = (hue + min_hue_separation) % 1.0

        # Add random perturbation to hue, saturation, and brightness
        hue = (hue + np.random.uniform(-0.02, 0.02)) % 1.0
        saturation = min(max(saturation + np.random.uniform(-0.05, 0.05), 0.0), 1.0)
        brightness = min(max(brightness + np.random.uniform(-0.05, 0.05), 0.0), 1.0)

        child_hsv = np.array([hue, saturation, brightness])
        child_color = to_hex(hsv_to_rgb(child_hsv))

        # Ensure uniqueness by perturbing hue slightly if the color is already used
        while child_color in used_colors:
            hue = (hue + min_hue_separation) % 1.0
            child_hsv[0] = hue
            child_color = to_hex(hsv_to_rgb(child_hsv))

        child_colors.append(child_color)
        used_colors.add(child_color)

    return child_colors


def _with_ancestors(marker_dict: dict, cell_types: set[str]) -> set[str]:
    """The cell types and their ancestors in the hierarchy."""
    out = set()

    def traverse(node, path):
        for cell_type, value in node.items():
            if cell_type in cell_types:
                out.update(path, {cell_type})
            traverse(value.get("subtypes", {}), [*path, cell_type])

    traverse(marker_dict, [])
    return out


def _build_hierarchy(
    marker_hierarchy: dict,
    present_cell_types: set[str],
    base_colormap: str,
    rotation_step: float,
    saturation_adjust: float,
    brightness_adjust: float,
    min_hue_separation: float,
) -> dict:
    """Restrict the marker hierarchy to the present cell types and give each a color derived from its parent's."""
    used_colors = set()

    def traverse(node, parent_color=None):
        hierarchy = {}
        children = [ct for ct in node.keys() if ct in present_cell_types]
        sibling_colors = (
            derive_child_colors(
                parent_color,
                len(children),
                rotation_step,
                saturation_adjust,
                brightness_adjust,
                min_hue_separation,
                used_colors,
            )
            if parent_color
            else []
        )

        for i, cell_type in enumerate(children):
            color = top_level_colors[cell_type] if parent_color is None else sibling_colors[i]
            subtypes = traverse(node[cell_type].get("subtypes", {}), parent_color=color)
            hierarchy[cell_type] = {"color": color, "subtypes": subtypes} if subtypes else {"color": color}

        return hierarchy

    top_level_colormap = mpl.colormaps.get_cmap(base_colormap)
    top_level_cell_types = [ct for ct in marker_hierarchy.keys() if ct in present_cell_types]
    num_top_level = len(top_level_cell_types)

    if hasattr(top_level_colormap, "colors"):
        available_colors = len(top_level_colormap.colors)
        if available_colors < num_top_level:
            raise ValueError(
                f"The colormap '{base_colormap}' only has {available_colors} colors, "
                f"but {num_top_level} are required for the top-level categories."
            )
        colors = [to_hex(color) for color in top_level_colormap.colors[:num_top_level]]
    else:
        # Sample continuous colormaps evenly
        colors = [to_hex(top_level_colormap(i / max(1, num_top_level - 1))) for i in range(num_top_level)]

    top_level_colors = dict(zip(top_level_cell_types, colors, strict=True))

    return traverse(marker_hierarchy)


def _assign_colors_and_order(hierarchy: dict) -> dict[str, str]:
    """Flatten the colored hierarchy depth-first into an ordered category-to-color mapping."""
    color_mapping = {}
    used_colors = set()

    def traverse(node):
        for cell_type, value in node.items():
            color = value["color"]
            if color in used_colors:
                raise ValueError(f"Duplicate color detected: {color} for cell type '{cell_type}'")
            used_colors.add(color)
            color_mapping[cell_type] = color
            if "subtypes" in value:
                traverse(value["subtypes"])

    traverse(hierarchy)
    return color_mapping


def process_snapseed(
    adata: ad.AnnData,
    marker_dict: dict,
    ann_df: pd.DataFrame,
    obs_key: str = "leiden",
    prefix: str = "snapseed",
    fill_na: bool = True,
    base_colormap: str = "tab10",
    get_colors: bool = True,
    sample_key: str | None = None,
    min_abs_freq: int | None = 5,
    min_rel_freq: float | None = 1e-3,
    rotation_step: float = 0.05,
    saturation_adjust: float = 0.2,
    brightness_adjust: float = 0.1,
    min_hue_separation: float = 0.05,
) -> None:
    """
    Map the cluster assignments of `annotate_hierarchy()` onto `adata.obs`, one column per level.

    Categories follow the hierarchy depth-first; colors are derived from the parent's color, so that
    subtypes look like their parent across levels.

    Args:
        adata: AnnData object.
        marker_dict: The marker hierarchy passed to `annotate_hierarchy()`.
        ann_df: The `"assignments"` returned by `annotate_hierarchy()`: one column per level
            (`level_1`, `level_2`, ...), indexed by cluster.
        obs_key: Column in `adata.obs` with the cluster labels passed to `annotate_hierarchy()`.
        prefix: Prefix for the new columns, named `<prefix>_level_<n>`.
        fill_na: Whether clusters without an annotation at a level keep that of the previous (coarser) level.
        base_colormap: Name of a matplotlib colormap for the top-level categories.
        get_colors: Whether to write colors consistent across levels to `adata.uns["<prefix>_level_<n>_colors"]`.
        sample_key: Column in `adata.obs` with the sample. If given, per sample and level, categories below
            `min_abs_freq` cells or below a fraction `min_rel_freq` of the sample's annotated cells are set to missing.
        min_abs_freq: Minimum number of cells of a category in a sample to retain it.
        min_rel_freq: Minimum fraction of a sample's annotated cells in a category to retain it.
        rotation_step: Step size for hue rotation when generating child colors.
        saturation_adjust: Adjustment factor for saturation (positive to increase, negative to decrease).
        brightness_adjust: Adjustment factor for brightness (positive to increase, negative to decrease).
        min_hue_separation: Minimum separation between hue values of child colors.

    Returns
    -------
        Nothing; writes `adata.obs["<prefix>_level_<n>"]` and, with `get_colors`, their colors to `adata.uns`.
    """
    levels = sorted(ann_df.keys(), key=lambda x: int(x.split("_")[1]))

    # a subtype can pass `min_expr` where its parent does not; the parent still needs a color to derive the subtype's
    present_cell_types = _with_ancestors(marker_dict, set(ann_df.stack().unique()))
    hierarchy = _build_hierarchy(
        marker_dict,
        present_cell_types,
        base_colormap=base_colormap,
        rotation_step=rotation_step,
        saturation_adjust=saturation_adjust,
        brightness_adjust=brightness_adjust,
        min_hue_separation=min_hue_separation,
    )
    cat_to_color = _assign_colors_and_order(hierarchy)
    global_order = list(cat_to_color.keys())

    if fill_na:
        ann_df = ann_df.ffill(axis=1)

    logger.info("Iterating over %s annotation levels to map to adata.obs", len(levels))
    for level in levels:
        col_name = f"{prefix}_{level}"
        adata.obs[col_name] = adata.obs[obs_key].map(ann_df[level].to_dict()).astype("category")

        if sample_key and (min_abs_freq or min_rel_freq):
            filtered = adata.obs[col_name].copy()
            n_cells_filtered = 0

            for sample, group in adata.obs.groupby(sample_key, observed=False):
                category_counts = group[col_name].value_counts()
                total_count = category_counts.sum()

                filtered_categories = set()
                if min_abs_freq:
                    filtered_categories.update(category_counts[category_counts < min_abs_freq].index)
                if min_rel_freq:
                    filtered_categories.update(category_counts[(category_counts / total_count) < min_rel_freq].index)

                sample_mask = adata.obs[sample_key] == sample
                cell_mask = ~filtered[sample_mask].isin(filtered_categories)
                filtered[sample_mask] = filtered[sample_mask].where(cell_mask, None)
                n_cells_filtered += np.sum(~cell_mask)

            logger.info(
                "Filtering: assigned %s/%s cells (%.2f%%) at '%s' with None.",
                n_cells_filtered,
                f"{adata.n_obs:,}",
                100 * (n_cells_filtered / adata.n_obs),
                level,
            )
            adata.obs[col_name] = filtered.cat.remove_unused_categories()

        # Order categories as in the hierarchy
        current_categories = adata.obs[col_name].cat.categories
        adata.obs[col_name] = adata.obs[col_name].cat.reorder_categories(
            [cat for cat in global_order if cat in current_categories]
        )

        if get_colors:
            adata.uns[f"{col_name}_colors"] = [cat_to_color[cat] for cat in adata.obs[col_name].cat.categories]
