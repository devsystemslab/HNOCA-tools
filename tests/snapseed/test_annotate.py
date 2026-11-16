import pandas as pd
import pytest

from hnoca.snapseed.annotate import annotate, annotate_hierarchy
from hnoca.snapseed.utils import get_annot_df, marker_dict_depth


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


def test_annotate(adata_annotate, marker_h1, expected_snap_class):
    """Test that annotate works with a non-hierarchical marker dict."""
    ann_df = annotate(adata_annotate, marker_h1, group_name="leiden")
    mapping = ann_df.set_index("leiden")["class"].to_dict()
    adata_annotate.obs["snap_class"] = adata_annotate.obs["leiden"].map(mapping)

    assert_annotation_counts(adata_annotate, "snap_class", expected_snap_class, "snap_class counts")


def test_annotate_no_markers_in_data(adata_annotate):
    """Marker genes that are all absent from `var_names` (e.g. Ensembl IDs vs. symbols) raise a clear error."""
    markers = {"T_cell": ["ENSG00000010610"], "B_cell": ["ENSG00000177455"]}
    with pytest.raises(ValueError, match="None of the marker genes"):
        annotate(adata_annotate, markers, group_name="leiden")


def test_annotate_some_markers_missing(adata_annotate):
    """Marker genes absent from `var_names` are ignored."""
    markers = {"T_cell": ["CD3E", "NOT_A_GENE"], "B_cell": ["MS4A1"]}
    assert set(annotate(adata_annotate, markers, group_name="leiden")["class"]) <= {"T_cell", "B_cell"}


def test_annotate_hierarchy(adata_annotate, marker_h2, expected_snap_level_1, expected_snap_level_2):
    """Test that annotate_hierarchy works with a hierarchical marker dict."""
    ann_dict = annotate_hierarchy(adata_annotate, marker_h2, group_name="leiden")
    ann_df = ann_dict["assignments"]

    # Map each level from the hierarchical annotation
    for key in ann_df:
        adata_annotate.obs[f"snap_{key}"] = adata_annotate.obs["leiden"].map(ann_df[key].to_dict())

    assert_annotation_counts(adata_annotate, "snap_level_1", expected_snap_level_1, "snap_level_1 counts")
    assert_annotation_counts(adata_annotate, "snap_level_2", expected_snap_level_2, "snap_level_2 counts")


def test_get_annot_df_call_needs_parent():
    """A call stands only below a call that stands, whatever its own detection rate."""

    def level(classes, expr):
        return pd.DataFrame({"leiden": ["0", "1", "2"], "class": classes, "expr": expr})

    metrics = {
        "level_1": level(["A", "A", "B"], [0.5, 0.01, 0.5]),
        "level_2": level(["A1", "A1", "B1"], [0.5, 0.5, 0.01]),
        "level_3": level(["A1a", "A1a", "B1a"], [0.5, 0.5, 0.5]),
    }
    ann = get_annot_df(metrics, "leiden", min_expr=0.1)
    assert ann.loc["0"].tolist() == ["A", "A1", "A1a"]
    # cluster 1 fails at level 1 and cluster 2 at level 2: nothing below stands, although it passes min_expr
    assert ann.loc["1"].isna().all()
    assert ann.loc["2", "level_1"] == "B"
    assert ann.loc["2", ["level_2", "level_3"]].isna().all()

    # without a threshold every call stands
    assert get_annot_df(metrics, "leiden", min_expr=0).notna().all().all()
