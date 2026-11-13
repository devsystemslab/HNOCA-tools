import anndata as ad
import numpy as np
import pandas as pd
import pytest

from hnoca.mapping.wknn import estimate_presence_score, get_wknn, transfer_labels


@pytest.fixture
def ref_query():
    rng = np.random.default_rng(0)
    ref = ad.AnnData(
        obs=pd.DataFrame(
            {"celltype": pd.Categorical(rng.choice(["a", "b", "c"], 60))}, index=[f"r{i}" for i in range(60)]
        ),
        obsm={"X_latent": rng.normal(size=(60, 5))},
    )
    query = ad.AnnData(
        obs=pd.DataFrame({"sample": pd.Categorical(rng.choice(["s1", "s2"], 40))}, index=[f"q{i}" for i in range(40)]),
        obsm={"X_latent": rng.normal(size=(40, 5))},
    )
    return ref, query


def test_transfer_labels(ref_query):
    ref, query = ref_query
    wknn = get_wknn(ref=ref.obsm["X_latent"], query=query.obsm["X_latent"], k=10)
    scores = transfer_labels(ref, query, wknn, label_key="celltype")
    assert list(scores.columns[-2:]) == ["best_label", "best_score"]
    assert scores["best_label"].isin(["a", "b", "c"]).all()
    np.testing.assert_allclose(scores["best_score"], scores[["a", "b", "c"]].max(axis=1))


def test_presence_score_split_by(ref_query):
    ref, query = ref_query
    res = estimate_presence_score(ref, query, k_wknn=10, split_by="sample", do_random_walk=False)
    assert set(res["per_group"].columns) == {"s1", "s2"}
    assert res["per_group"].index.equals(ref.obs_names)
