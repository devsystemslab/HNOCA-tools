from importlib.resources import files

import pandas as pd
import pytest
import scanpy as sc
import yaml


@pytest.fixture
def marker_h1():
    data_path = files("tests.data") / "snapseed_markers_h1.yaml"
    return yaml.safe_load(data_path.read_text())


@pytest.fixture
def marker_h2():
    data_path = files("tests.data") / "snapseed_markers_h2.yaml"
    return yaml.safe_load(data_path.read_text())


@pytest.fixture
def marker_h3():
    data_path = files("tests.data") / "snapseed_markers_h3.yaml"
    return yaml.safe_load(data_path.read_text())


@pytest.fixture
def adata_annotate():
    adata = sc.datasets.pbmc3k()

    # basic cell and gene filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Saving count data
    adata.layers["counts"] = adata.X.copy()

    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)

    # compute hvgs
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # PCA
    sc.tl.pca(adata, mask_var="highly_variable")

    # k-NN and leiden clustering
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution=5, random_state=0, directed=False)

    return adata


# Ground-truth Series fixtures
@pytest.fixture
def expected_snap_class():
    return pd.Series([348, 713, 141, 1498], index=["B_cell", "Monocyte", "NK_cell", "T_cell"], name="snap_class")


@pytest.fixture
def expected_snap_level_1():
    return pd.Series([348, 713, 141, 1498], index=["B_cell", "Monocyte", "NK_cell", "T_cell"], name="snap_level_1")


@pytest.fixture
def expected_snap_level_2():
    return pd.Series(
        [1322, 176, 542, 171],
        index=["CD4_T_cell", "CD8_T_cell", "Classical_monocyte", "Non_classical_monocyte"],
        name="snap_level_2",
    )
