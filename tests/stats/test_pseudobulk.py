import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from hnoca.stats.pseudobulk import create_pseudobulk


@pytest.fixture
def adata_counts():
    rng = np.random.default_rng(0)
    n = 120
    X = sparse.csr_matrix(rng.poisson(5, size=(n, 30)).astype(np.float32))
    obs = pd.DataFrame(
        {
            "sample": pd.Categorical(np.repeat(["s1", "s2", "s3"], n // 3)),
            "celltype": pd.Categorical(np.tile(["a", "b"], n // 2)),
        },
        index=[f"c{i}" for i in range(n)],
    )
    return ad.AnnData(X=X, obs=obs, layers={"counts": X.copy()})


@pytest.mark.parametrize("group_key", [None, "celltype"])
def test_create_pseudobulk(adata_counts, group_key):
    obs_before = adata_counts.obs.copy()
    pb = create_pseudobulk(adata_counts, sample_key="sample", group_key=group_key, layer="counts", min_counts=10)

    assert pb.n_obs == (3 if group_key is None else 6)
    np.testing.assert_allclose(pb.X.sum(), adata_counts.X.sum())
    expected = pd.Series(np.asarray((adata_counts.X > 0).sum(axis=1)).ravel()).groupby(obs_before["sample"].values)
    np.testing.assert_allclose(pb.obs["n_genes_median"], pb.obs["sample"].map(expected.median()).astype(float))
    # the input is not modified
    pd.testing.assert_frame_equal(adata_counts.obs, obs_before)


def test_create_pseudobulk_filters(adata_counts):
    pb = create_pseudobulk(adata_counts, sample_key="sample", group_key="celltype", min_cells=100)
    assert pb.n_obs == 0


def test_create_pseudobulk_no_filter(adata_counts):
    pb = create_pseudobulk(adata_counts, sample_key="sample", group_key="celltype", min_cells=0, min_counts=0)
    assert pb.n_obs == 6
