import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import ListedColormap, rgb_to_hsv, to_hex, to_rgb

from hnoca.snapseed import annotate_hierarchy, display_missing_genes, process_snapseed
from hnoca.snapseed.process import derive_child_colors

TAB10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]


@pytest.fixture
def ann_df():
    """Assignments as returned by annotate_hierarchy() for the depth-3 PBMC markers, with missing calls."""
    return pd.DataFrame(
        {
            "level_1": ["T_cell", "T_cell", "T_cell", "B_cell", "Monocyte", np.nan],
            "level_2": ["CD4_T_cell", "CD4_T_cell", "CD8_T_cell", np.nan, "Classical_monocyte", np.nan],
            "level_3": ["Treg", np.nan, np.nan, np.nan, np.nan, np.nan],
        },
        index=pd.Index([str(i) for i in range(6)], name="leiden"),
    )


@pytest.fixture
def adata_clusters():
    """Two samples: cluster 0 has 3 cells in sample A and 1 in sample B, all other clusters 10 cells per sample."""
    sizes = {"A": [3, 10, 10, 10, 10, 10], "B": [1, 10, 10, 10, 10, 10]}
    obs = pd.DataFrame(
        [(str(c), s) for s, ns in sizes.items() for c, n in enumerate(ns) for _ in range(n)],
        columns=["leiden", "sample"],
    )
    obs.index = obs.index.astype(str)
    return ad.AnnData(obs=obs.astype("category"))


def test_process_snapseed_fill_na(adata_clusters, marker_h3, ann_df):
    process_snapseed(adata_clusters, marker_h3, ann_df, prefix="snap", min_abs_freq=None, min_rel_freq=None)
    obs = adata_clusters.obs
    first = obs.groupby("leiden", observed=True)[["snap_level_1", "snap_level_2", "snap_level_3"]].first()

    # cluster 1 keeps its level-2 call at level 3, cluster 3 its level-1 call at levels 2 and 3
    assert first.loc["1"].tolist() == ["T_cell", "CD4_T_cell", "CD4_T_cell"]
    assert first.loc["3"].tolist() == ["B_cell", "B_cell", "B_cell"]
    assert first.loc["0", "snap_level_3"] == "Treg"
    # a cluster without any call stays missing
    assert obs.loc[obs["leiden"] == "5", "snap_level_3"].isna().all()


def test_process_snapseed_no_fill(adata_clusters, marker_h3, ann_df):
    process_snapseed(adata_clusters, marker_h3, ann_df, prefix="snap", fill_na=False, sample_key=None)
    obs = adata_clusters.obs
    assert obs.loc[obs["leiden"] == "3", "snap_level_2"].isna().all()
    assert obs["snap_level_3"].cat.categories.tolist() == ["Treg"]


def test_process_snapseed_order_and_colors(adata_clusters, marker_h3, ann_df):
    np.random.seed(0)
    process_snapseed(adata_clusters, marker_h3, ann_df, prefix="snap", sample_key=None)
    uns, obs = adata_clusters.uns, adata_clusters.obs

    # categories follow the hierarchy depth-first, absent cell types are left out
    assert obs["snap_level_1"].cat.categories.tolist() == ["T_cell", "B_cell", "Monocyte"]
    assert obs["snap_level_2"].cat.categories.tolist() == [
        "CD4_T_cell",
        "CD8_T_cell",
        "B_cell",
        "Classical_monocyte",
    ]

    # top-level colors come from the colormap in hierarchy order, and a cell type has the same color at every level
    assert uns["snap_level_1_colors"] == TAB10[:3]
    colors = {}
    for level in ["snap_level_1", "snap_level_2", "snap_level_3"]:
        assert len(uns[f"{level}_colors"]) == len(obs[level].cat.categories)
        for cat, color in zip(obs[level].cat.categories, uns[f"{level}_colors"], strict=True):
            assert colors.setdefault(cat, color) == color
    assert len(set(colors.values())) == len(colors)


def test_process_snapseed_sample_filter(adata_clusters, marker_h3, ann_df):
    # sample B has a single Treg cell (below min_abs_freq), sample A has 3 (kept)
    process_snapseed(adata_clusters, marker_h3, ann_df, prefix="snap", sample_key="sample", min_abs_freq=2)
    obs = adata_clusters.obs
    treg = obs["leiden"] == "0"
    assert obs.loc[treg & (obs["sample"] == "A"), "snap_level_3"].eq("Treg").all()
    assert obs.loc[treg & (obs["sample"] == "B"), "snap_level_3"].isna().all()

    # min_rel_freq: 3 of 53 cells in sample A fall below a fraction of 0.1
    process_snapseed(
        adata_clusters, marker_h3, ann_df, prefix="rel", sample_key="sample", min_abs_freq=None, min_rel_freq=0.1
    )
    assert "Treg" not in obs["rel_level_3"].cat.categories
    assert len(adata_clusters.uns["rel_level_3_colors"]) == len(obs["rel_level_3"].cat.categories)


def test_process_snapseed_annotate_hierarchy(adata_annotate, marker_h3):
    ann_df = annotate_hierarchy(adata_annotate, marker_h3, group_name="leiden")["assignments"]
    process_snapseed(adata_annotate, marker_h3, ann_df, prefix="snap")

    for level in ann_df:
        assert f"snap_{level}" in adata_annotate.obs
        assert f"snap_{level}_colors" in adata_annotate.uns
    # with fill_na, deeper levels only refine: every cell's level-2 label is its level-1 label or a subtype of it
    children = {k: set(v.get("subtypes", {})) | {k} for k, v in marker_h3.items()}
    obs = adata_annotate.obs.dropna(subset=["snap_level_1", "snap_level_2"])
    assert all(l2 in children[l1] for l1, l2 in zip(obs["snap_level_1"], obs["snap_level_2"], strict=True))


def test_process_snapseed_child_without_parent(adata_clusters, marker_h3, ann_df):
    # snapseed recurses before filtering by min_expr, so a cluster can pass at level 2 but not at level 1
    ann_df.loc["5", "level_2"] = "CD8_T_cell"
    process_snapseed(adata_clusters, marker_h3, ann_df, prefix="snap", sample_key=None)
    obs = adata_clusters.obs
    cl5 = obs["leiden"] == "5"
    assert obs.loc[cl5, "snap_level_1"].isna().all()
    assert obs.loc[cl5, "snap_level_2"].eq("CD8_T_cell").all()
    assert len(adata_clusters.uns["snap_level_2_colors"]) == len(obs["snap_level_2"].cat.categories)


def test_process_snapseed_small_colormap(adata_clusters, marker_h3, ann_df):
    # three top-level cell types, two colors
    name = "hnoca_test_two_colors"
    if name not in mpl.colormaps:
        mpl.colormaps.register(ListedColormap(["#000000", "#ffffff"], name=name))
    with pytest.raises(ValueError, match="only has 2 colors, but 3 are required"):
        process_snapseed(adata_clusters, marker_h3, ann_df, base_colormap=name, sample_key=None)


def test_derive_child_colors():
    np.random.seed(0)
    used = {"#123456"}
    children = derive_child_colors("#1f77b4", 4, 0.05, 0.2, 0.1, 0.05, used_colors=used)
    assert len(set(children)) == 4
    assert used == {"#123456", *children}
    parent_hue = rgb_to_hsv(np.array(to_rgb("#1f77b4")))[0]
    assert all(to_hex(c) != "#1f77b4" for c in children)
    assert all(abs(rgb_to_hsv(np.array(to_rgb(c)))[0] - parent_hue) > 0.02 for c in children)

    assert derive_child_colors("#1f77b4", 0, 0.05, 0.2, 0.1, 0.05) == []
    with pytest.raises(ValueError, match="Parent color cannot be None"):
        derive_child_colors(None, 2, 0.05, 0.2, 0.1, 0.05)


def test_display_missing_genes(monkeypatch, adata_annotate, marker_h3):
    messages = []
    monkeypatch.setattr("hnoca.snapseed.process.logger.info", lambda msg, *args: messages.append(msg % args))

    display_missing_genes(marker_h3, adata_annotate)
    missing = messages[1:]
    assert messages[0] == "Genes not found in adata.var:"
    assert "TBET" in missing
    assert missing == sorted(missing)
    assert not set(missing) & set(adata_annotate.var_names)

    messages.clear()
    display_missing_genes({"B_cell": {"marker_genes": ["MS4A1", "CD79A"]}}, adata_annotate)
    assert messages == ["All marker genes are present in adata.var."]
