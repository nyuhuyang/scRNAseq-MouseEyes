# scRNAseq-MouseEyes

## Research Background
The choroid is a highly vascularized layer of the eye localized between the sclera and the outermost retinal layer, the retinal pigment epithelium (RPE). Blood supplied by choroidal circulation is the main source of oxygen and nutrients for the RPE and photoreceptors, as well as the main evacuation route for retinal waste. On the other hand, the RPE provides essential support functions for photoreceptor homeostasis and visual function. Thus, it is not surprising that RPE/choroid alterations are present in a vast range of posterior ophthalmic diseases including uveal inflammation and macular degenerative diseases such as Best disease, Stargardt disease and age-related macular degeneration (AMD).

AMD affects 25% people over 80 years old, being a leading cause of blindness in developed countries and approaching epidemic proportions in the United States. However, currently, there is the only treatment for 10% of patients who suffer the advanced, angiogenic form of the disease, and not all respond successfully. The lack of therapeutic options is mainly due to the fact that AMD etiology remains unknown, and probably involves the malfunction of several cell types and multiple intercellular signaling pathways. Thus, it is imperative to carry out systematic studies to characterize in detail all RPE/choroid cell types at the molecular level and to understand how cells communicate among them to maintain choroid and retinal homeostasis. This would constitute the first step to assess whether choroidal intercellular cellular crosstalk is compromised in AMD and to develop new therapeutic strategies based on the restoration of such signaling circuits.

Here, we used adult mice to carry out for the first time single cell RNAseq of RPE/choroid tissue, providing its molecular characterization at single cell resolution. In combination with tissue-specific endothelial cell transcription profiling and the use of transgenic mouse models, we discovered a novel immunomodulatory signaling circuit in the choroid.

## Data
Chromium single-cell RNA-seq outputs were processed by Cell Ranger analysis pipelines. Data is currently unavailable to the public before publication.

## How to use this repository

### Software Setup
R version 3.4.3 (Did't test other versions)<br />
dplyr_0.7.4 (Did't test other versions)<br />
Seurat_2.2.1 (must be >2.2.0 )<br />

After pull this repository, create folders **_data_** and **_output_** in current working folder.
Move Cell Ranger analysis results in to **_data_** folder.

### Seurat_setup.R
Unsupervised cell clustering analysis was carried out using the Seurat 2.2 R package. Cells with <500 genes and genes detected within <3 cells were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and a log transformation, using the Seurat NormalizeData function. The top 1000 highly variable genes from young C57BL/6J and aged C57BL/6J datasets were selected, followed by a canonical correlation analysis (CCA) to identify common sources of variation between the two datasets and minimize the batch effect. The first 20 CCA results were chosen for principal component analysis (PCA). Cells were used for 2-dimensional t-Distributed Stochastic Neighbor Embedding (tSNE) (ref van der maaten and hinton 2008) with 0.8 resolution.

 After running this script, the `mouse_eyes_alignment.Rda` file will be generated inside **_data_** folder.
 Do not modify any files in **_data_** folder.
 
 
### Identify_Cell_Types_Manually.R
All clusters are examed against 122(number may change) CD marker genes.
All cell types are predicted by at least two marker genes with adjusted p-value smaller than 10^-30.

Endothelial cells were identified by Cdh5, Flt1, Kdr, Pecam1, Plvap, Ptprb, and Vwf.<br />
Pericytes were identified by Dcn, Des, Ifitm1, Mylk, Pdgfrb, and Rgs5.<br />
Hematopoietic cells were identified by Laptm5, Ptprc, and Srgn.<br />
Melanocytes were identified by Mlana and Pmel.<br />
Myelinating Schwann cells were identified by Mbp and Mpz.<br />
Retinal pigment epitheliums were identified by Rlbp1 and Rpe65.<br />


Cells contained in cluster 11 (hematopoietic cells) were further subjected to a second round of unsupervised analysis following the same approach, resulting in a tSNE analysis with ~0.1 resolution. The modified Seurat function FindAllMarkers was used to calculate average differential expression among cell clusters. The p-value was calculated using likelihood-ratio test and adjusted by Benjamini-Hochberg method.
