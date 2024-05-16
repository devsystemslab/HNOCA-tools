import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .tests import ancova_group_test, f_nonzero_test, anova


def test_de(
    adata,
    group,
    covar=None,
    num_threads=1,
    return_coef_group=None,
    var_names=None,
    adjust_method="holm",
):
    """
    Test for differential expression using ANOVA

    Parameters
    ----------
    adata : AnnData
        The AnnData object
    group : str or pd.Series
        The group labels
    covar : str or pd.DataFrame
        The covariates
    num_threads : int
        The number of threads to use
    return_coef_group : str
        The group to return the coefficient for
    var_names : list
        The variable names to test
    adjust_method : str
        The method to adjust p-values. See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

    Returns
    -------
    pd.DataFrame
        The differential expression results
    """
    if var_names is None:
        var_names = adata.var_names

    expr_mat = adata[:, var_names].X

    if isinstance(covar, str):
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
    query_adata,
    matched_adata,
    covar=None,
    num_threads=1,
    var_names=None,
    adjust_method="holm",
):
    """
    Test for differential expression between query data and matches reference cells using an F-test.

    Parameters
    ----------
    query_adata : AnnData
        The query AnnData object
    matched_adata : AnnData
        The matched reference AnnData object
    covar : str or pd.DataFrame
        The covariates
    num_threads : int
        The number of threads to use
    var_names : list
        The variable names to test
    adjust_method : str
        The method to adjust p-values. See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    """
    if var_names is None:
        var_names = query_adata.var_names
    else:
        var_names = np.intersect1d(var_names, query_adata.var_names)

    var_names = np.intersect1d(var_names, matched_adata.var_names)

    if isinstance(covar, str):
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
