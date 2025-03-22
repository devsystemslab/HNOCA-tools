import pandas as pd
import pytest

from hnoca.snapseed.annotate import annotate, annotate_hierarchy
from hnoca.snapseed.utils import marker_dict_depth


def assert_annotation_counts(adata, obs_column, expected_series, obj_str):
    """Helper to assert that the counts in adata.obs[obs_column] match expected_series."""
    actual = adata.obs.groupby(obs_column, observed=True).size()
    actual.index = actual.index.astype(object)
    pd.testing.assert_series_equal(
        actual.sort_index(),
        expected_series.sort_index(),
        obj=obj_str,
        check_index_type=False,
        check_names=False,
    )


@pytest.mark.parametrize(
    "marker_fixture, expected_depth",
    [
        ("marker_h1", 1),
        ("marker_h2", 2),
        ("marker_h3", 3),
    ],
)
def test_marker_depth(request, marker_fixture, expected_depth):
    # Retrieve the fixture value by name using request.getfixturevalue
    markers = request.getfixturevalue(marker_fixture)
    depth = marker_dict_depth(markers)
    assert depth == expected_depth, f"Expected depth {expected_depth} for fixture {marker_fixture}, got {depth}"


def test_annotate_hierarchy_invalid_depth(adata_annotate, marker_h1):
    """Test that annotate_hierarchy raises a ValueError when passed a non-hierarchical marker dict."""
    with pytest.raises(ValueError, match="You passed a non-hierachical marker dict to annotate_hierarchy"):
        annotate_hierarchy(adata_annotate, marker_h1, group_name="leiden")


@pytest.mark.parametrize("marker_fixture", ["marker_h2", "marker_h3"])
def test_annotate_invalid_depth(request, adata_annotate, marker_fixture):
    """Test that annotate raises a ValueError when passed a hierarchical marker dict."""
    markers = request.getfixturevalue(marker_fixture)
    with pytest.raises(ValueError, match="Use annotate_hierarchy"):
        annotate(adata_annotate, markers, group_name="leiden")


def test_annotate(adata_annotate, marker_h1):
    """Test that annotate works with a non-hierarchical marker dict."""
    ann_df = annotate(adata_annotate, marker_h1, group_name="leiden")
    mapping = ann_df.set_index("leiden")["class"].to_dict()
    adata_annotate.obs["snap_class"] = adata_annotate.obs["leiden"].map(mapping)

    expected_counts = pd.Series(
        [348, 713, 141, 1498], index=["B_cell", "Monocyte", "NK_cell", "T_cell"], name="snap_class"
    )

    assert_annotation_counts(adata_annotate, "snap_class", expected_counts, "snap_class counts")


def test_annotate_hierarchy(adata_annotate, marker_h2):
    """Test that annotate_hierarchy works with a hierarchical marker dict."""
    ann_dict = annotate_hierarchy(adata_annotate, marker_h2, group_name="leiden")
    ann_df = ann_dict["assignments"]

    # Map each level from the hierarchical annotation
    for key in ann_df:
        adata_annotate.obs[f"snap_{key}"] = adata_annotate.obs["leiden"].map(ann_df[key].to_dict())

    expected_level1 = pd.Series(
        [348, 713, 141, 1498], index=["B_cell", "Monocyte", "NK_cell", "T_cell"], name="snap_level_1"
    )
    assert_annotation_counts(adata_annotate, "snap_level_1", expected_level1, "snap_level_1 counts")

    expected_level2 = pd.Series(
        [1322, 176, 542, 171],
        index=["CD4_T_cell", "CD8_T_cell", "Classical_monocyte", "Non_classical_monocyte"],
        name="snap_level_2",
    )
    assert_annotation_counts(adata_annotate, "snap_level_2", expected_level2, "snap_level_2 counts")
