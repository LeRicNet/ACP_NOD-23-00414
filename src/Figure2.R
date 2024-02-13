# Figure 2 Code
library(Seurat)
library(tidyverse)
library(reshape2)

spatial.obj <- readRDS('./data/seurat_objects/ACP_visium.rds')
spatial.obj$Br

# Figure 2A
SpatialFeaturePlot(
  spatial.obj,
  features = c(
    'Brown', ## Apps Epithelial Module
    'Magenta', ## Apps Immune Module
    'Blue' ## Apps Glial Module
  ),
  images = c('slice1', 'slice1_sampleB')
)

# Figure 2B
FeatureScatter(
  spatial.obj,
  'Brown',
  'Blue',
  pt.size = 4
) +
  labs(
    y = 'glial',
    x = 'epithelialCTNNB1 mutation',
    title = ''
  ) +
  theme_linedraw(base_size=14)

## Figure 2C
# This figure is generated using results from Metascape.org. The gene lists that
# are used as inputs for metascape can be derived with the following code:
epi.glial.markers <- FindMarkers(spatial.obj,
                           ident.1 = 'Epithelial',
                           ident.2 = 'NonImmune')
epi.markers <- epi.glial.markers %>%
  rownames_to_column('gene') %>%
  as_tibble() %>%
  filter(avg_log2FC > 0,
         p_val_adj < 1e-3) %>%
  arrange(desc(abs(avg_log2FC)))

glial.markers <- epi.glial.markers %>%
  rownames_to_column('gene') %>%
  as_tibble() %>%
  filter(avg_log2FC < 0,
         p_val_adj < 1e-3) %>%
  arrange(desc(abs(avg_log2FC)))


## Figure 2D
spatial.obj@meta.data %>%
  as_tibble() %>%
  select(sample, Brown, Blue, Magenta) %>%
  melt(id.vars=c('sample')) %>%
  ggplot(aes(variable, value, fill=variable)) +
  geom_boxplot() +
  facet_wrap(~sample, nrow = 1) +
  theme_bw(base_size=14)
