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
