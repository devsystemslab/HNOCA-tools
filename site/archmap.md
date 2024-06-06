# How to map query data to the HNOCA using ArchMap

Here we give a short overview of the relevant features of the [ArchMap](https://www.archmap.bio) web service.

Archmap is an easy-to-use web frontend that offers query-to-reference mapping for scRNA-seq datasets and interactive exploration of the results. It is, therefore, a great resource for mapping your new organoid dataset to the HNOCA or the fetal brain atlas by [Braun et al.](https://doi.org/10.1126/science.adf1226)

Please refer to the official ArchMap documentation for a general introduction and explanation of how to use ArchMap for query-to-reference mapping.


## Use case 1: Mapping an organoid query dataset to the HNOCA with ArchMap

For quick contextualisation and annotation of your organoid datasets, you can use ArchMap to conveniently map it to the HNOCA. After creating your account, simply select the HNOCA from the ArchMap core reference atlases, select the scPoli model and a classifier (we recommend the KNN classifier for quick results), and upload your dataset. Please note the requirements listed on the submit page and adapt your adata object accordingly before uploading.

Once the upload and processing has been completed, you can either download the mapped dataset or interactively explore the results through the inbuilt CellxGene interface.

The following obs fields in the results object might be of interest:

- *predictions_knn*: here you can find the transferred HNOCA cell type labels for your query dataset
- *type*: this annotation lets you quickly distinguish your query dataset from the reference data
- *Hallmark_Glycolysis*: in the HNOCA paper, we introduce the concept of a metabolic stress score, which is specific to organoid datasets. ArchMap automatically computes this stress score for your query cells
- *\*_uncertainty_euclidean*: this provides an uncertainty score per cell regarding the cell type label prediction found in *predictions_knn*
- *presence_score*: as introduced in the HNOCA paper, this score indicates, for every **reference** cell, how well that cell state is represented in your query dataset.


## Use case 2: Mapping an organoid query dataset to the fetal brain atlas

If you want to contextualise your organoid dataset using the recent fetal brain atlas by [Braun et al.](https://doi.org/10.1126/science.adf1226), you can follow the same steps as described above but select the fetal brain reference instead of HNOCA from the ArchMap Core Atlases.

**Note:** Mapping organoid data to a primary reference is a challenging integration task, and the mapping might not be satisfactory depending on the query data. We are working on improving this.

The following obs fields in the results object might be of interest:

- *predictions_knn*: here you can find the transferred HNOCA cell type labels for your query dataset
- *type*: this annotation lets you quickly distinguish your query dataset from the reference data
- *\*_uncertainty_euclidean*: this provides an uncertainty score per cell regarding the cell type label prediction found in *predictions_knn*
- *presence_score*: as introduced in the HNOCA paper, this score indicates, for every **reference** cell, how well that cell state is represented in your query dataset.
