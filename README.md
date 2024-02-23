# ACP_NOD-23-00414

## Overview

This repository contains computational analyses and models on Adamantinomatous Craniopharyngioma (ACP), a serious CNS tumor known for significantly affecting patient's quality of life. Current treatment only extends to surgery and radiation. Our research centers on investigating the role of senescence and the Senescence-Associated Secretory Phenotype (SASP) in ACP primary tumor tissue. For additional information, please see:

> Prince EW, Apps JR, Jeang J, Chee K, Medlin S, Jackson EM, Dudley R, Limbrick D, Naftel R, Johnston J, Feldstein N, Prolo LM, Ginn K, Niazi T, Smith A, Kilburn L, Chern J, Leonard J, Lam S, Hersh DS, Gonzalez-Meljem JM, Amani V, Donson AM, Mitra SS, Bandohpadhayay P, Martinez-Barbera JP, Hankinson TC. **Unraveling the Complexity of the Senescence-Associated Secretory Phenotype in Adamantinomatous Craniopharyngioma Using Multi-Modal Machine Learning Analysis.** Neuro Oncol. 2024 Feb 9:noae015. doi: 10.1093/neuonc/noae015. <a href='https://pubmed.ncbi.nlm.nih.gov/38334125/'>PMID: 38334125</a>.

## Resources

This repository is divided into two components, Data and Code. The data in this repository includes the Seurat objects for the single-nucleus, -cell, and spatial gene expression data sets. The code in this repository is organized by Figure and Figure Panel in the manuscript text.

<hr/>

### Code to Generate Paper Figures

<p align="center"><a href="src/Figure1.R">Figure 1</a> |
 <a href="src/Figure2.R">Figure 2</a> |
 <a href="src/Figure3.R">Figure 3</a> |
 <a href="src/Figure4.R">Figure 4</a> |
 <a href="src/Figure5.R">Figure 5</a> |
 <a href="src/Figure6.R">Figure 6</a>
</p>

<hr/>

### Data

For more information see [here](docs/appendix.md).

`data/seurat_objects/acp_scn_annotated.rds` (\~1.3 GB)

`data/geneset_enrichment_matrices/acp_scn_msigdb_em.rds` (\~1.7 GB)

`data/seurat_objects/acp_spatial_annotated.rds` (\~353 MB)




## Troubleshooting

Please report any issues, comments, or questions to Eric Prince via email [Eric.Prince\@CUAnschutz.edu](mailto:Eric.Prince@CUAnschutz.edu), or [file an issue](https://github.com/LeRicNet/ACP_NOD-23-00414/issues).
