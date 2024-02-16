# ACP_NOD-23-00414

## Overview

This repository contains computational analyses and models on Adamantinomatous Craniopharyngioma (ACP), a serious CNS tumor known for significantly affecting patient's quality of life. Current treatment only extends to surgery and radiation. Our research centers on investigating the role of senescence and the Senescence-Associated Secretory Phenotype (SASP) in ACP primary tumor tissue. For additional information, please see:

> Prince EW, Apps JR, Jeang J, Chee K, Medlin S, Jackson EM, Dudley R, Limbrick D, Naftel R, Johnston J, Feldstein N, Prolo LM, Ginn K, Niazi T, Smith A, Kilburn L, Chern J, Leonard J, Lam S, Hersh DS, Gonzalez-Meljem JM, Amani V, Donson AM, Mitra SS, Bandohpadhayay P, Martinez-Barbera JP, Hankinson TC. **Unraveling the Complexity of the Senescence-Associated Secretory Phenotype in Adamantinomatous Craniopharyngioma Using Multi-Modal Machine Learning Analysis.** Neuro Oncol. 2024 Feb 9:noae015. doi: 10.1093/neuonc/noae015. Epub ahead of print. <a href='https://pubmed.ncbi.nlm.nih.gov/38334125/'>PMID: 38334125</a>.

## Resources

This repository is divided into two components, Data and Code. The data in this repository includes the Seurat objects for the single-nucleus, -cell, and spatial gene expression data sets. The code in this repository is organized by Figure and Figure Panel in the manuscript text.


### Data

`data/seurat_objects/ACP_singlecell-nucleus.rds` (~1.6 GB)

`data/seurat_objects/ACP_visium.rds` (~355 MB)


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
