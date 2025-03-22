import pandas as pd
import pytest

from hnoca.snapseed.annotate import annotate, annotate_hierarchy
from hnoca.snapseed.utils import marker_dict_depth


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
    """
    Test that annotate_hierarchy raises a ValueError when passed a non-hierarchical marker dict.
    """
    with pytest.raises(ValueError, match="You passed a non-hierachical marker dict to annotate_hierarchy"):
        annotate_hierarchy(adata_annotate, marker_h1, group_name="leiden")


@pytest.mark.parametrize("marker_fixture", ["marker_h2", "marker_h3"])
def test_annotate_invalid_depth(request, adata_annotate, marker_fixture):
    """
    Test that annotate raises a ValueError when passed a hierarchical marker dict.
    """
    markers = request.getfixturevalue(marker_fixture)
    with pytest.raises(ValueError, match="Use annotate_hierarchy"):
        annotate(adata_annotate, markers, group_name="leiden")


def test_annotate_h1(adata_annotate, marker_h1):
    """
    Test that annotate works with a non-hierarchical marker dict.
    """
    # Generate annotations with snapseed
    ann_df = annotate(adata_annotate, marker_h1, group_name="leiden")
    mapping = ann_df.set_index("leiden")["class"].to_dict()
    adata_annotate.obs["snapseed_class"] = adata_annotate.obs["leiden"].map(mapping)

    # Get the actual counts
    actual_counts = adata_annotate.obs.groupby("snapseed_class", observed=True).size()
    actual_counts.index = actual_counts.index.astype(object)

    # Create the expected pandas Series
    expected_counts = pd.Series(
        [348, 713, 141, 1498], index=["B_cell", "Monocyte", "NK_cell", "T_cell"], name="snapseed_class"
    )

    # Compare the to Series
    pd.testing.assert_series_equal(
        actual_counts.sort_index(),
        expected_counts.sort_index(),
        obj="snapseed_class counts",
        check_index_type=False,
        check_names=False,
    )
