# ACP_NOD-23-00414

## Overview

This repository contains computational analyses and models on Adamantinomatous Craniopharyngioma (ACP), a serious CNS tumor known for significantly affecting patient's quality of life. Current treatment only extends to surgery and radiation. Our research centers on investigating the role of senescence and the Senescence-Associated Secretory Phenotype (SASP) in ACP primary tumor tissue. For additional information, please see:

> Prince EW, Apps JR, Jeang J, Chee K, Medlin S, Jackson EM, Dudley R, Limbrick D, Naftel R, Johnston J, Feldstein N, Prolo LM, Ginn K, Niazi T, Smith A, Kilburn L, Chern J, Leonard J, Lam S, Hersh DS, Gonzalez-Meljem JM, Amani V, Donson AM, Mitra SS, Bandohpadhayay P, Martinez-Barbera JP, Hankinson TC. **Unraveling the Complexity of the Senescence-Associated Secretory Phenotype in Adamantinomatous Craniopharyngioma Using Multi-Modal Machine Learning Analysis.** Neuro Oncol. 2024 Feb 9:noae015. doi: 10.1093/neuonc/noae015. <a href='https://pubmed.ncbi.nlm.nih.gov/38334125/'>PMID: 38334125</a>.

## Resources

This repository is divided into two components, Data and Code. The data in this repository includes the Seurat objects for the single-nucleus, -cell, and spatial gene expression data sets. The code in this repository is organized by Figure and Figure Panel in the manuscript text.

### Data

`data/seurat_objects/acp_scn_annotated.rds` (\~1.3 GB)

This RDS file contains a [Seurat](https://satijalab.org/seurat/) object that stores single-cell RNA-seq data from the ACP project. The data has been preprocessed, normalized, and integrated using the [Seurat pipeline](https://satijalab.org/seurat/archive/v3.2/integration). The object has the following characteristics:

-   It has 48,814 features (genes) across 18,900 samples (10,322 cells & 8,578 nuclei) within three assays: integrated, RNA, and SCT.
-   The active assay is integrated (cells and nuclei), which contains the data after integration and dimensionality reduction. It has 3,000 features, of which 2,195 are variable features that capture the heterogeneity of the data. Data has been integrated according to the Seurat protocol.
-   The integrated assay has two layers: data and scale.data. The data layer contains the log-normalized expression values and the scale.data layer includes the scaled and centered expression values.
-   The other two assays are RNA and SCT, which contain the raw and corrected expression values. They have the same number of features and samples as the integrated assay but different layers.
-   The object has two dimensional reductions: PCA and UMAP. The PCA reduction contains the integrated data's principal components (PCs), which are used for clustering and finding marker genes. The UMAP reduction contains the uniform manifold approximation and projection (UMAP) coordinates of the integrated data, which are used for visualization and exploration.

The metadata for this object has the following columns:

-   `level`: The level of sequencing used for the Seurat object. It can be either "cell" or "nucleus," indicating whether the data was obtained from single cells or single nuclei.
-   `nCount_RNA`: The number of unique molecular identifiers (UMIs) assigned to each cell or nucleus in the RNA assay. It reflects the sequencing depth and quality of the sample.
-   `nFeature_RNA`: The number of features (genes) detected in each cell or nucleus in the RNA assay. It reflects the complexity and diversity of the transcriptome of the sample.
-   `sample_id`: The identifier of the sample that the cell or nucleus belongs to. It can be used to group or subset the data by sample origin.
-   `nCount_SCT`: The number of UMIs assigned to each cell or nucleus in the SCT assay. It is similar to nCount_RNA, but after applying the SCTransform normalization method.
-   `nFeature_SCT`: The number of features detected in each cell or nucleus in the SCT assay. It is similar to nFeature_RNA, but after applying the SCTransform normalization method.
-   `seurat_clusters`: The Seurat pipeline assigns the cluster labels to each cell or nucleus. They are based on the integrated data's PCA and UMAP dimensional reductions and reflect the samples' similarities and dissimilarities.
-   `celltypes`: The cell type annotations assigned to each cell or nucleus by the SignacX algorithm. They are based on the gene expression pattern of the samples, and map them to the Human Cell Landscape, a comprehensive atlas of human cell types.


`data/genesets/msigdb.rds`

This file contains a `GeneSetCollection` object that stores gene sets from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb). The gene sets are typically used for gene set enrichment analysis (GSEA) or other signature-based analyses. The object has the following characteristics:

- It has 25,377 gene sets from the MSigDB, covering various biological processes, pathways, and phenotypes. The gene sets are named according to the MSigDB conventions, such as HALLMARK_ADIPOGENESIS, HALLMARK_ALLOGRAFT_REJECTION, etc.
- It has 40,525 unique identifiers for the genes in the gene sets. The identifiers are gene symbols, such as ABCA1 and ABCB8.
- The object was created using the [escape](https://www.bioconductor.org/packages/devel/bioc/manuals/escape/man/escape.pdf) package from Bioconductor, which provides an interface to access the MSigDB in R. The R call to get the gene set is: `msigdb <- getGeneSets(species = "Homo sapiens", library = c("H", paste0("C", 1:6)))`. This call retrieves the gene sets for the human species and the libraries H (hallmark gene sets) and C1-C6 (curated gene sets).

`data/geneset_enrichment_matrices/acp_scn_msigdb_em.rds` (\~1.7 GB)

This file contains an R `data.frame` object with 25,377 gene sets (columns) and 18,900 samples (rows). The values are the enrichment scores for each of the samples (i.e., cells/nuclei) and gene sets. The object was created using the `escape` package, using the call `ES <- enrichIt(obj = s.obj.raw, gene.sets = msigdb, method = "UCell", cores = 32, maxRank = 10000)`.

`data/seurat_objects/acp_spatial_annotated.rds` (\~353 MB)

This Seurat object contains single-cell transcriptomic and spatial data from four samples of the 10X Visium platform. The samples are named slice1, slice1_sampleB, slice1_sampleC, and slice1_sampleD, and their corresponding images are stored in the object. The object has three assays: SCT, Spatial, and RNA. The SCT assay contains the normalized and scaled data for 20,322 features (genes), of which 3000 are variable across 8,366 cells (spots). The Spatial and RNA assays contain the same raw count data and the spatial coordinates for each spot. The object also has four dimensional reductions: pca, umap, harmony, and umap_harmony, which can be used to cluster and visualize the cells.

The metadata for this object has the following columns:

-   `sample_id`: The identifier of the sample that the spot belongs to. It can be used to group or subset the data by sample origin.
-   `nCount_Spatial`: The number of unique molecular identifiers (UMIs) assigned to each spot in the RNA assay. It reflects the sequencing depth and quality of the sample.
-   `nFeature_Spatial`: The number of features (genes) detected in each spot in the RNA assay. It reflects the complexity and diversity of the transcriptome of the sample.
-   `nCount_SCT`: The number of UMIs assigned to each spot in the SCT assay. It is similar to nCount_RNA, but after applying the SCTransform normalization method.
-   `nFeature_SCT`: The number of features detected in each spot in the SCT assay. It is similar to nFeature_RNA, but after applying the SCTransform normalization method.
-   `seurat_clusters`: The Seurat pipeline assigns the cluster labels to each spot or tissue. They are based on the integrated data's PCA and UMAP dimensional reductions and reflect the samples' similarities and dissimilarities.
-   `cellstates`: The cell type annotations to each spot by the SignacX algorithm. They are based on the gene expression pattern of the samples and map them to the Human Cell Landscape, a comprehensive atlas of human cell types.


### Code

`src/Figure1.R`

`src/Figure2.R`

`src/Figure3.R`

`src/Figure4.R`

`src/Figure5.R`

`src/Figure6.R`

### Exporting Figures and Postprocessing

All figure panels were exported as EPS files and formatted using Adobe Illustrator. To make code more readable, the exportation of figure panels has been omitted. The process to export a figure panel, using Figure 2A as an example, would be:

```         
setEPS()

postscript(
  file = "./figure_2a.eps",
  height = 7,
  width = 4
)

# Figure 2A (found in src/Figure2.R)
SpatialFeaturePlot(
  spatial.obj,
  features = c(
    'Brown', ## Apps Epithelial Module
    'Magenta', ## Apps Immune Module
    'Blue' ## Apps Glial Module
  ),
  images = c('slice1', 'slice1_sampleB')
)

dev.off()
```

## Troubleshooting

Please report any issues, comments, or questions to Eric Prince via email [Eric.Prince\@CUAnschutz.edu](mailto:Eric.Prince@CUAnschutz.edu), or [file an issue](https://github.com/LeRicNet/ACP_NOD-23-00414/issues).
