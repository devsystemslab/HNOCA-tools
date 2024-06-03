from typing import Union, Optional
import anndata as ad
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .tests import ancova_group_test, f_nonzero_test, anova


def test_de(
    adata: ad.AnnData,
    group: Union[str, pd.Series],
    covar: Union[str, pd.DataFrame],
    num_threads: int = 1,
    return_coef_group: str = None,
    var_names: list = None,
    adjust_method: str = "holm",
) -> pd.DataFrame:
    """
    Test for differential expression using ANOVA

    Args:
        adata: AnnData object
        group: str or pd.Series
            The group labels
        covar: str or pd.DataFrame
            The covariates
        num_threads: int
            The number of threads to use
        return_coef_group: str
            The group to return coefficients for
        var_names: list
            The variable names to test
        adjust_method: str
            The method to adjust p-values. See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

    Returns:
        A `pd.DataFrame` with the differential expression results
    """
    if var_names is None:
        var_names = adata.var_names

    expr_mat = adata[:, var_names].X

    if isinstance(covar, (str, list, tuple)):
        covar = adata.obs[covar]
    elif isinstance(covar, pd.Series):
        covar = covar.values

    if isinstance(group, str):
        group = adata.obs[group]
    elif isinstance(group, pd.Series):
        group = group.values

    results = ancova_group_test(
        expr_mat,
        group,
        covar=covar,
        num_threads=num_threads,
        return_coef_group=return_coef_group,
        var_names=var_names,
    )

    results["padj"] = multipletests(results["p_Resi"], method=adjust_method)[1]
    results["pval"] = results["p_Resi"]
    return results


def test_de_paired(
    query_adata: ad.AnnData,
    matched_adata: ad.AnnData,
    covar: Optional[Union[str, pd.DataFrame]] = None,
    num_threads: int = 1,
    var_names: list = None,
    adjust_method: str = "holm",
) -> pd.DataFrame:
    """
    Test for differential expression between query data and matches reference cells using an F-test.

    Args:
        query_adata: AnnData object
            The query data
        matched_adata: AnnData object
            The matched reference data
        covar: str or pd.DataFrame
            The covariates
        num_threads: int
            The number of threads to use
        var_names: list
            The variable names to test
        adjust_method: str
            The method to adjust p-values. See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

    Returns:
        A `pd.DataFrame` with the differential expression results
    """
    if var_names is None:
        var_names = query_adata.var_names
    else:
        var_names = np.intersect1d(var_names, query_adata.var_names)

    var_names = np.intersect1d(var_names, matched_adata.var_names)

    if isinstance(covar, (str, list, tuple)):
        covar = query_adata.obs[covar]
    elif isinstance(covar, pd.Series):
        covar = covar.values

    expr_0 = query_adata[:, var_names].X - matched_adata[:, var_names].X

    results = f_nonzero_test(
        expr_0,
        covar=covar,
        num_threads=num_threads,
        var_names=var_names,
    )

    results["padj"] = multipletests(results["pval"], method=adjust_method)[1]
    return results
