##**************************************************************************##
##*
##* Figure 3: Transcriptional signature of SASP cells in ACP epithelial tissue.
##*
##**************************************************************************##

# Load Dependencies
library(Seurat)
library(tidyverse)
library(escape)

# Load single-cell, -nucleus gene expression data
s.obj <- readRDS("./data/seurat_objects/acp_scn_annotated.rds")

# Load spatial gene expression data
spatial.obj <- readRDS("./data/seurat_objects/acp_spatial_annotated.rds")

##**************************************************************************##
##* Figure 3A
##**************************************************************************##

# Subset epithelial samples
s.obj <- subset(s.obj, celltypes == 'Epithelial')

# Calculate enrichment for GenAge, CellAge, and SenMayo genesets
age_gene_sets <- readRDS("./data/genesets/aging-and-senescence.rds")
ES <- enrichIt(
  obj = s.obj,
  gene.sets = age_gene_sets
)
s.obj <- AddMetaData(s.obj, scale(ES))

# Label SASP(+)/(-) based on the 90th quantile threshold
s.obj@meta.data$SASP <- ifelse(s.obj$SenMayo > quantile(s.obj$SenMayo, 0.9),
                               'SASP(+)', 'SASP(-)')
s.obj@meta.data %>%
  ggplot(aes(CellAge, GenAge, col=SASP)) +
  geom_point()

##**************************************************************************##
##* Figure 3B
##**************************************************************************##


##**************************************************************************##
##* Figure 3C
##**************************************************************************##
##*
##* Calculate the enrichment scores for the ageing gene sets within the
##* saptial gene expression data.

ES2 <- enrichIt(
  obj = spatial.obj,
  gene.sets = age_gene_sets
)
spatial.obj <- AddMetaData(spatial.obj, scale(ES2))
spatial.obj@meta.data$SASP <- ifelse(
  spatial.obj$SenMayo > quantile(spatial.obj$SenMayo, 0.9),
  'SASP(+)', 'SASP(-)')

Idents(spatial.obj) <- 'SASP'
SpatialDimPlot(spatial.obj, images = 'slice1')

##**************************************************************************##
##* Figure 3D
##**************************************************************************##
##*
##* Extract the differentially expressed genes (DEGs) within the epithelial cells
##* and nuclei that correspond to the top and bottom 10% of the SenMayo score.
##* These genes were then exported to Metascape.org to determine which gene
##* sets, pathways, and ontologies were associated with a given gene set.

s.obj@meta.data$label <- ifelse(s.obj$SenMayo > quantile(s.obj$SenMayo, 0.9),
                                'top10', NA)
s.obj@meta.data$label <- ifelse(s.obj$SenMayo < quantile(s.obj$SenMayo, 0.1),
                                'bottom10', s.obj$label)
Idents(s.obj) <- 'label'
top10_senmayo_degs <- FindMarkers(s.obj, ident.1 = 'top10') # should return
                                                            # n=1,307 genes/rows
bottom10_senmayo_degs <- FindMarkers(s.obj, ident.1 = 'bottom10')

##**************************************************************************##
##* Figure 3E
##**************************************************************************##
##* TODO

##**************************************************************************##
##* Figure 3F
##**************************************************************************##
##* TODO
