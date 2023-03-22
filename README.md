# Cross-species single-cell comparison of systemic and cardiac inflammatory responses after cardiac injury

Cortada et al., 2023. <https://www.biorxiv.org/content/10.1101/2023.03.15.532865v1.full>

>Scripts for the analysis of scRNA-seq data of mouse and zebrafish 

## Methods

Various tissues from both zebrafish and mouse were prepared for scRNA-seq to identify the organism-wide gene signatures following experimental myocardial infarct (MI) in fish and mice.

### scRNA-seq analysis

The raw reads were aligned and processed with the CellRanger pipeline (v6.0.0) using the mouse (mm10) or zebrafish genomes (GRCz11).\
Subsequent analyses were performed in R following the recommendations of Amezquita et al. (<https://osca.bioconductor.org/>) using numerous functions provided in the R packages  `scater` and `scran`. 
Briefly, quality control was carried out for each sample separately with functions from the `scuttle` package; cells with low gene content and high mitochondrial gene content were removed from further analyses.
For each species, different count matrices across all samples were scaled to account for sequencing depth differences and log-transformed using `multiBatchNorm` from the `batchelor` package. Cell types were annotated with `SingleR` using numerous datasets and manual inspection of marker genes.
Data was also explored and annotated with `Cellxgene` version v0.16.8. 
scRNA-seq data UMAPs, SAMaps, dot plots, and violin plots were generated with Cellxgene VIP (<https://doi.org/10.1101/2020.08.28.270652>).
<
#### Cross-species comparisons

`SAMap` (v1.0.2) was used to determine cell homology between mouse and zebrafish myeloid clusters.
Briefly, a reciprocal BLAST map between mouse and zebrafish transcriptomes was generated using the `map_genes.sh` script from SAMap, which ran `tblastx` (nucleotide-nucleotide BLAST, NCBI) in both directions with respect to Ensembl zebrafish and mouse transcriptomes (`Danio_rerio.GRCz11.cds.all.fa` and `Mus_musculus.GRCm38.cds.all.fa`).
Raw scRNA-seq count matrices of zebrafish and mouse were extracted from `SingleCellExperiment` objects with cluster annotations and converted to `h5ad` format using the `sceasy` R package (v0.0.7). Raw data were then processed and integrated by `SAMap`. Mapping scores between the myeloid clusters of each species were calculated with SAMap `get_mapping_scores` (n top=0) and visualized in R using `ggalluvial`; cluster pairs with alignment scores less than 0.2 were excluded.

For myeloid cluster pairs with high inter-species alignment scores (mm1-dr4, mm4-dr6, and mm6-dr3), marker genes for each cluster were assessed using scran::scoreMarkers (mean AUC > 0.6), and over-representation analyses were performed using `clusterProfiler` with respect to gene ontology terms. To compare the transcriptional response to MI (mouse) or cardiac cryoinjury (zebrafish) between the myeloid cluster pairs, differential expression analysis was performed using scran::findMarkers(lfc=0.25) at days 1 and 7 in the heart and in the blood. Genes were considered statistically significant at an FDR cutoff of 0.10.

## References

* scater: McCarthy, D.J., Campbell, K.R., Lun, A.T. & Wills, Q.F. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics 33, 1179–1186 (2017).
* Lun, A.T., McCarthy, D.J. & Marioni, J.C. A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Res 5, 2122 (2016)
* SingleR: Aran, D. et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat Immunol 20, 163–172 (2019).
* MNN: Haghverdi, L., Lun, A.T.L., Morgan, M.D. & Marioni, J.C. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat Biotechnol 36, 421–427 (2018).
* clusterProfiler: Yu, G., Wang, L.G., Han, Y. & He, Q.Y. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16, 284–287 (2012).
* SAMap: Tarashansky, A.J. et al. Mapping single-cell atlases throughout Metazoa unravels cell type evolution. Elife 10 (2021)
