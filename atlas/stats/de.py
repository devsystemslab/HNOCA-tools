import pandas as pd
from statsmodels.stats.multitest import multipletests

from .tests import ancova_group_test, f_nonzero_test, anova


def test_de(
    adata, group, covar=None, num_threads=1, return_coef_group=None, var_names=None
):
    """
    Test for differential expression using ANOVA

    Parameters
    ----------
    adata : AnnData
        The AnnData object
    group : pd.Series
        The group labels
    covar : pd.DataFrame
        The covariates
    num_threads : int
        The number of threads to use
    return_coef_group : str
        The group to return the coefficient for
    var_names : list
        The variable names to test

    Returns
    -------
    pd.DataFrame
        The differential expression results
    """
    if var_names is None:
        var_names = adata.var_names

    expr_mat = adata[:, var_names].X

    if covar is not None:
        covar = adata.obs[covar]

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

    results["padj"] = multipletests(results["p_Resi"], method="holm")[1]
