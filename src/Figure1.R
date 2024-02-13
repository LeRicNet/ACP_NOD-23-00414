# Figure 1 Code
library(Seurat)
library(tidyverse)

# Load Visium spatial gene expression Seurat Object
spatial.obj <- readRDS('./data/seurat_objects/ACP_visium.rds')

# Set default identities as CellStates from Human Cell Landscape
Idents(spatial.obj) <- 'CellStates'

# Figure 1A contains Visium slide images

# Figure 1B
DimPlot(spatial.obj, reduction='umap_harmony', pt.size = 4)

# Figure 1C
SpatialDimPlot(spatial.obj)

# Figure 1D
Idents(spatial.obj) <- 'seurat_clusters'
DimPlot(spatial.obj, reduction='umap_harmony', pt.size = 4)

# Figure 1E
SpatialDimPlot(spatial.obj)

# Figure 1F; TODO: Revise this entry; there is a discrepency between original
#            figure and this one. Might be a color ramp issue.
spatial.obj@meta.data %>%
  as_tibble() %>%
  group_by(sample, seurat_clusters, CellStates) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(seurat_clusters, CellStates, fill=pct)) +
  geom_tile() +
  facet_wrap(~sample, nrow=1)

# Figure 1G
# This figure is generated using results from Metascape.org. The gene lists that
# are used as inputs for metascape can be derived with the following code:

# TODO: this is also not running..
spatially.variable.features <- FindSpatiallyVariableFeatures(spatial.obj,
                                                             assay = 'RNA',
                                                             slot = 'counts',
                                                             selection.method = 'moransi')
