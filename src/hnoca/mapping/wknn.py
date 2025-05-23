import warnings
from typing import Literal

import numpy as np
import pandas as pd
import scanpy as sc
import tqdm
from pynndescent import NNDescent
from scipy import sparse

from hnoca._logging import logger
from hnoca.utils.check import check_deps


def nn2adj_gpu(nn, n1=None, n2=None):  # noqa: D103
    if n1 is None:
        n1 = nn[1].shape[0]
    if n2 is None:
        n2 = np.max(nn[1].flatten())

    df = pd.DataFrame(
        {
            "i": np.repeat(range(nn[1].shape[0]), nn[1].shape[1]),
            "j": nn[1].flatten(),
            "x": nn[0].flatten(),
        }
    )
    adj = sparse.csr_matrix((np.repeat(1, df.shape[0]), (df["i"], df["j"])), shape=(n1, n2))

    return adj


def nn2adj_cpu(nn, n1=None, n2=None):  # noqa: D103
    if n1 is None:
        n1 = nn[0].shape[0]
    if n2 is None:
        n2 = np.max(nn[0].flatten())

    df = pd.DataFrame(
        {
            "i": np.repeat(range(nn[0].shape[0]), nn[0].shape[1]),
            "j": nn[0].flatten(),
            "x": nn[1].flatten(),
        }
    )

    adj = sparse.csr_matrix((np.repeat(1, df.shape[0]), (df["i"], df["j"])), shape=(n1, n2))

    return adj


def build_nn(  # noqa: D103
    ref,
    query=None,
    k=100,
    use_rapids: bool = False,
):
    if query is None:
        query = ref

    if use_rapids:
        check_deps("cuml")
        from cuml.neighbors import NearestNeighbors

        logger.info("Using cuML for neighborhood estimation on GPU.")

        model = NearestNeighbors(n_neighbors=k)
        model.fit(ref)
        knn = model.kneighbors(query)
        adj = nn2adj_gpu(knn, n1=query.shape[0], n2=ref.shape[0])
    else:
        logger.info("Using pynndescent for neighborhood estimation on CPU.")
        index = NNDescent(ref)
        knn = index.query(query, k=k)
        adj = nn2adj_cpu(knn, n1=query.shape[0], n2=ref.shape[0])

    return adj


def build_mutual_nn(dat1, dat2=None, k1=100, k2=None):  # noqa: D103
    if dat2 is None:
        dat2 = dat1
    if k2 is None:
        k2 = k1

    adj_12 = build_nn(dat1, dat2, k=k2)
    adj_21 = build_nn(dat2, dat1, k=k1)

    adj_mnn = adj_12.multiply(adj_21.T)
    return adj_mnn


def get_transition_prob_mat(dat, k=50, symm=True):  # noqa: D103
    adj = build_nn(dat, k=k)
    if symm:
        adj = adj + adj.transpose()
    prob = sparse.diags(1 / np.array(adj.sum(1)).flatten()) @ adj.transpose()
    return prob


def random_walk_with_restart(init, transition_prob, alpha=0.5, num_rounds=100):  # noqa: D103
    init = np.array(init).flatten()
    heat = init[:, None]
    for _i in range(num_rounds):
        heat = init[:, None] * alpha + (1 - alpha) * (transition_prob.transpose() @ heat)
    return heat


def get_wknn(
    ref: np.ndarray,
    query: np.ndarray,
    ref2: np.ndarray | None = None,
    k: int = 100,
    query2ref: bool = True,
    ref2query: bool = False,
    weighting_scheme: Literal["n", "top_n", "jaccard", "jaccard_square", "gaussian", "dist"] = "jaccard_square",
    top_n: int | None = None,
    return_adjs: bool = False,
):
    """
    Compute the weighted k-nearest neighbors graph between the reference and query datasets

    Parameters
    ----------
    ref
        The reference representation to build ref-query neighbor graph
    query
        The query representation to build ref-query neighbor graph
    ref2
        The reference representation to build ref-ref neighbor graph
    k
        Number of neighbors per cell
    query2ref
        Consider query-to-ref neighbors
    ref2query
        Consider ref-to-query neighbors
    weighting_scheme
        How to weight edges in the ref-query neighbor graph
    top_n
        The number of top neighbors to consider
    return_adjs
        Whether to return the adjacency matrices of ref-query, query-ref, ref-ref, and ref-ref for weighting
    """
    adj_q2r = build_nn(ref=ref, query=query, k=k)

    adj_r2q = None
    if ref2query:
        adj_r2q = build_nn(ref=query, query=ref, k=k)

    if query2ref and not ref2query:
        adj_knn = adj_q2r.T
    elif ref2query and not query2ref:
        adj_knn = adj_r2q
    elif ref2query and query2ref:
        adj_knn = ((adj_r2q + adj_q2r.T) > 0) + 0
    else:
        warnings.warn(
            "At least one of query2ref and ref2query should be True. Reset to default with both being True.",
            stacklevel=2,
        )
        adj_knn = ((adj_r2q + adj_q2r.T) > 0) + 0

    if ref2 is None:
        ref2 = ref
    adj_ref = build_nn(ref=ref2, k=k)
    logger.info("Info: Computing shared neighbors")
    num_shared_neighbors = adj_q2r @ adj_ref.T
    logger.info("Info: Computing shared neighbors")
    num_shared_neighbors_nn = num_shared_neighbors.multiply(adj_knn.T)

    wknn = num_shared_neighbors_nn.copy()
    if weighting_scheme == "top_n":
        if top_n is None:
            top_n = k // 4 if k > 4 else 1
        wknn = (wknn > top_n) * 1
    elif weighting_scheme == "jaccard":
        wknn.data = wknn.data / (k + k - wknn.data)
    elif weighting_scheme == "jaccard_square":
        wknn.data = (wknn.data / (k + k - wknn.data)) ** 2

    if return_adjs:
        adjs = {"q2r": adj_q2r, "r2q": adj_r2q, "knn": adj_knn, "r2r": adj_ref}
        return (wknn, adjs)
    else:
        return wknn


def estimate_presence_score(
    ref_adata: sc.AnnData,
    query_adata: sc.AnnData,
    wknn=None,
    use_rep_ref_wknn: str = "X_latent",
    use_rep_query_wknn: str = "X_latent",
    k_wknn: int = 100,
    query2ref_wknn: bool = True,
    ref2query_wknn: bool = False,
    weighting_scheme_wknn: str = "jaccard_square",
    ref_trans_prop=None,
    use_rep_ref_trans_prop=None,
    k_ref_trans_prop: int = 50,
    split_by=None,
    do_random_walk: bool = True,
    alpha_random_walk: float = 0.1,
    num_rounds_random_walk: int = 100,
    log=True,
):
    """Estimate presence score of query cells in reference dataset."""
    if wknn is None:
        ref = ref_adata.obsm[use_rep_ref_wknn]
        query = query_adata.obsm[use_rep_query_wknn]
        wknn = get_wknn(
            ref=ref,
            query=query,
            k=k_wknn,
            query2ref=query2ref_wknn,
            ref2query=ref2query_wknn,
            weighting_scheme=weighting_scheme_wknn,
        )

    if ref_trans_prop is None and do_random_walk:
        if use_rep_ref_trans_prop is None:
            use_rep_ref_trans_prop = use_rep_ref_wknn
        ref = ref_adata.obsm[use_rep_ref_trans_prop]
        ref_trans_prop = get_transition_prob_mat(ref, k=k_ref_trans_prop)

    if split_by and split_by in query_adata.obs.columns:
        presence_split = [
            np.array(wknn[query_adata.obs[split_by] == x, :].sum(axis=0)).flatten()
            for x in query_adata.obs[split_by].unique()
        ]
    else:
        presence_split = [np.array(wknn.sum(axis=0)).flatten()]

    if do_random_walk:
        logger.info("Info: Smoothing presence scores with random walk")
        presence_split_sm = []
        for x in tqdm.tqdm(presence_split):
            presence_split_sm.append(
                random_walk_with_restart(
                    init=x,
                    transition_prob=ref_trans_prop,
                    alpha=alpha_random_walk,
                    num_rounds=num_rounds_random_walk,
                )
            )
    else:
        presence_split_sm = [x[:, None] for x in presence_split]

    columns = query_adata.obs[split_by].unique() if split_by and split_by in query_adata.obs.columns else ["query"]
    if len(columns) > 1:
        df_presence = pd.DataFrame(
            np.concatenate(presence_split_sm, axis=1),
            columns=columns,
            index=ref_adata.obs_names,
        )
    else:
        df_presence = pd.DataFrame({columns[0]: presence_split_sm[0].flatten()}).set_index(ref_adata.obs_names)

    if log:
        df_presence = df_presence.apply(lambda x: np.log1p(x), axis=0)
    df_presence_norm = df_presence.apply(lambda x: np.clip(x, np.percentile(x, 1), np.percentile(x, 99)), axis=0).apply(
        lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)), axis=0
    )
    max_presence = df_presence_norm.max(1)

    return {
        "max": max_presence,
        "per_group": df_presence_norm,
        "ref_trans_prop": ref_trans_prop,
    }


def transfer_labels(ref_adata: sc.AnnData, query_adata: sc.AnnData, wknn, label_key: str = "celltype"):
    """Transfer labels from reference to query data."""
    scores = pd.DataFrame(
        wknn @ pd.get_dummies(ref_adata.obs[label_key]),
        columns=pd.get_dummies(ref_adata.obs[label_key]).columns,
        index=query_adata.obs_names,
    )
    scores["best_label"] = scores.idxmax(1)
    scores["best_score"] = scores.max(1)
    return scores
