import pytest

from hnoca.snapseed.utils import marker_dict_depth


@pytest.mark.parametrize(
    "marker_fixture, expected_depth",
    [
        ("marker_hierarchy", 2),
        ("marker_flat", 1),
    ],
)
def test_marker_depth(request, marker_fixture, expected_depth):
    # Retrieve the fixture value by name using request.getfixturevalue
    markers = request.getfixturevalue(marker_fixture)
    depth = marker_dict_depth(markers)
    assert depth == expected_depth, f"Expected depth {expected_depth} for fixture {marker_fixture}, got {depth}"
