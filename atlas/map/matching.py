import anndata
import numpy as np
import scanpy as sc
from scipy import sparse


def get_matched_transcriptome(
    adata,
    adata_ref,
    wknn,
    rescale_factor=1,
):
    expr_ref = adata_ref.X.copy()
    if rescale_factor != 1:
        expr_ref = ((expr_ref.expm1()) * rescale_factor).log1p()

    normmat = sparse.diags(1 / np.array(wknn.sum(axis=1)).flatten())
    softmax_wknn = normmat.dot(wknn)
    expr_bg = softmax_wknn @ expr_ref

    adata_bg = anndata.AnnData(
        expr_bg, obs=adata.obs.copy(), var=adata_ref.var.copy(), obsm=adata.obsm.copy()
    )
    return adata_bg
