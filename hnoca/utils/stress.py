import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc


def compute_glycolysis_score(adata: ad.AnnData, **kwargs):
    """Compute glycolysis score using hallmark glycolysis genes from MSigDB.

    Args:
        adata: Anndata object with gene expression data.
        **kwargs: Additional arguments passed to `sc.tl.score_genes`.

    Returns:
        `None`. The glycolysis score is stored in `adata.obs["Hallmark_Glycolysis"]`.
    """
    msigdb_glycolysis = np.array(
        pd.read_csv(
            "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=HALLMARK_GLYCOLYSIS&fileType=TSV",
            sep="\t",
            header=None,
            index_col=0,
        )
        .loc["GENE_SYMBOLS", 1]
        .split(",")
    )
    msigdb_glycolysis = np.intersect1d(msigdb_glycolysis, adata.var_names)
    sc.tl.score_genes(
        adata, msigdb_glycolysis, score_name="Hallmark_Glycolysis", **kwargs
    )
