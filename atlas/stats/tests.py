import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import f
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from joblib import Parallel, delayed, wrap_non_picklable_objects
import multiprocessing
from tqdm import tqdm


def anova(e, group, covar=None, return_coef_group=None):
    if issparse(e):
        e = e.todense()
    e = np.array(e).flatten()

    if covar is None:
        m0 = ols("e ~ 1", data=pd.DataFrame({"e": e})).fit()
        m1 = ols("e ~ zz_group", data=pd.DataFrame({"e": e, "zz_group": group})).fit()
    else:
        data = pd.DataFrame({"e": e, **covar, "zz_group": group})
        m0 = ols("e ~ " + " + ".join(covar.columns), data=data).fit()
        m1 = ols("e ~ " + " + ".join(covar.columns) + " + zz_group", data=data).fit()

    a0 = anova_lm(m0)
    a1 = anova_lm(m1)
    a01 = anova_lm(m0, m1)
    p1 = a01.loc[:, "Pr(>F)"].iloc[-1]
    p2_cdf = f.cdf(
        a0.loc["Residual", "mean_sq"] / a1.loc["Residual", "mean_sq"],
        a0.loc["Residual", "df"],
        a1.loc["Residual", "df"],
    )
    p2 = 1 - np.abs(0.5 - p2_cdf) * 2

    var_tot = a0["sum_sq"].sum()
    var_covar = a0.iloc[:-1]["sum_sq"].sum()
    var_g = a1.loc["zz_group", "sum_sq"]

    if (
        return_coef_group is not None
        and return_coef_group in pd.Series(group).unique()
        and f"zz_group[T.{return_coef_group}]" in m1.params.index
    ):
        coef_g = m1.params[f"zz_group[T.{return_coef_group}]"]
    else:
        coef_g = np.nan

    return (var_g, var_covar, var_tot, p1, p2, coef_g)


delayed_anova = delayed(wrap_non_picklable_objects(anova))


def ancova_group_test(
    expr, group, covar=None, num_threads=1, return_coef_group=None, var_names=None
):
    if not isinstance(group, pd.Series):
        group = pd.Series(group)
    idx = np.where(~group.isnull())[0]
    expr = expr[idx, :]
    group = group[idx]
    if covar is not None:
        covar = pd.DataFrame(covar).iloc[idx, :]

    if group.dtype.name != "category":
        group = group.astype("category")
    group = group.cat.remove_unused_categories()
    g_lev = np.array(pd.Series(group).astype("category").cat.categories)
    if return_coef_group is not None and return_coef_group in g_lev:
        g_lev = np.concatenate([g_lev[g_lev != return_coef_group], [return_coef_group]])
        group = group.cat.reorder_categories(g_lev)

    num_cores = min(num_threads, multiprocessing.cpu_count())
    if num_cores > 1:
        results = Parallel(n_jobs=num_cores)(
            delayed_anova(expr[:, i], group, covar, return_coef_group)
            for i in range(expr.shape[1])
        )
    else:
        results = list()
        for i in tqdm(range(expr.shape[1])):
            results.append(
                anova(
                    expr[:, i], group, covar=covar, return_coef_group=return_coef_group
                )
            )

    df_res = pd.DataFrame(
        results,
        columns=["var_group", "var_covar", "var_total", "p_ANCOVA", "p_Resi", "coef"],
    )
    if var_names is not None and len(var_names) == df_res.shape[0]:
        df_res.index = var_names
    else:
        df_res.index = pd.RangeIndex(start=0, stop=expr.shape[1], step=1)
    return df_res


def f_nonzero(e, covar=None):
    if issparse(e):
        e = e.todense()
    e = np.array(e).flatten()

    if covar is None:
        m1 = ols("e ~ 1", data=pd.DataFrame({"e": e})).fit()
        a1 = anova_lm(m1)
        fval = (np.sum(e**2) / len(e)) / a1.loc["Residual", "mean_sq"]
        p = 1 - np.abs(0.5 - f.cdf(fval, len(e), a1.loc["Residual", "df"])) * 2
        coef = m1.params["Intercept"]
    else:
        data = pd.DataFrame({"e": e, **covar})
        m0 = ols("e ~ " + " + ".join(covar.columns) + " -1", data=data).fit()
        a0 = anova_lm(m0)
        m1 = ols("e ~ " + " + ".join(covar.columns), data=data).fit()
        a1 = anova_lm(m1)
        sumsq0 = a0.loc["Residual", "sum_sq"] + (
            a0.iloc[:-1]["sum_sq"].sum() - a1.iloc[:-1]["sum_sq"].sum()
        )
        meansq0 = sumsq0 / (a0.loc["Residual", "df"] + 1)
        meansq1 = a1.loc["Residual", "mean_sq"]
        fval = meansq0 / meansq1
        p = (
            1
            - np.abs(
                0.5
                - f.cdf(fval, a0.loc["Residual", "df"] + 1, a1.loc["Residual", "df"])
            )
            * 2
        )
        coef = m1.params["Intercept"]

    return (fval, coef, p)


delayed_f_nonzero = delayed(wrap_non_picklable_objects(f_nonzero))


def f_nonzero_test(expr, covar=None, num_threads=1, var_names=None):
    if covar is not None:
        covar = pd.DataFrame(covar)

    num_cores = min(num_threads, multiprocessing.cpu_count())
    if num_cores > 1:
        results = Parallel(n_jobs=num_cores)(
            delayed_f_nonzero(expr[:, i], covar) for i in range(expr.shape[1])
        )
    else:
        results = list()
        for i in tqdm(range(expr.shape[1])):
            results.append(f_nonzero(expr[:, i], covar=covar))

    df_res = pd.DataFrame(results, columns=["f", "coef", "pval"])
    if var_names is not None and len(var_names) == df_res.shape[0]:
        df_res.index = var_names
    else:
        df_res.index = pd.RangeIndex(start=0, stop=expr.shape[1], step=1)
    return df_res
