# scRNAseq-MouseEyes

 | ----- | ----------------------------------------
 | Title | Single cell transcriptomics of adult and young mouse RPE/choroid 
 | Author | Yang Hu
 | Department | Physiology and Biophysics
 | University | Weill Cornell Medicine
 | Email | yah2014@med.cornell.edu


## Research Background
The choroid is a highly vascularized layer of the eye localized between the sclera and the outermost retinal layer, the retinal pigment epithelium (RPE). Blood supplied by choroidal circulation is the main source of oxygen and nutrients for the RPE and photoreceptors, as well as the main evacuation route for retinal waste. On the other hand, the RPE provides essential support functions for photoreceptor homeostasis and visual function. Thus, it is not surprising that RPE/choroid alterations are present in a vast range of posterior ophthalmic diseases including uveal inflammation and macular degenerative diseases such as Best disease, Stargardt disease and age-related macular degeneration (AMD).

AMD affects 25% people over 80 years old, being a leading cause of blindness in developed countries and approaching epidemic proportions in the United States. However, currently, there is the only treatment for 10% of patients who suffer the advanced, angiogenic form of the disease, and not all respond successfully. The lack of therapeutic options is mainly due to the fact that AMD etiology remains unknown, and probably involves the malfunction of several cell types and multiple intercellular signaling pathways. Thus, it is imperative to carry out systematic studies to characterize in detail all RPE/choroid cell types at the molecular level and to understand how cells communicate among them to maintain choroid and retinal homeostasis. This would constitute the first step to assess whether choroidal intercellular cellular crosstalk is compromised in AMD and to develop new therapeutic strategies based on the restoration of such signaling circuits.

Here, we used young and aged mice to carry out for the first time single cell RNAseq of RPE/choroid tissue, providing its molecular characterization at single cell resolution. In combination with tissue-specific endothelial cell transcription profiling and the use of transgenic mouse models, we discovered a novel immunomodulatory signaling circuit in the choroid.

## Data
Chromium single-cell RNA-seq outputs were processed by Cell Ranger analysis pipelines. Data is currently unavailable to the public before publication.

## How to use this repository

#### Software Setup
R version 3.4.3 (Did't test other versions)<br />
dplyr_0.7.4 (Did't test other versions)<br />
Seurat_2.2.1 (Must be >2.2.0 )<br />

After pulling this repository, create folders **_data_** and **_output_** in the top working folder.
Move Cell Ranger analysis results into **_data_** folder.

### 1. Seurat_setup.R
Unsupervised cell clustering analysis was carried out using the Seurat 2.2 R package. Cells with <500 genes and genes detected within <3 cells were excluded from the analysis. Gene expression raw counts were normalized following a global-scaling normalization method with a scale factor of 10,000 and a log transformation, using the Seurat NormalizeData function. The top 1000 highly variable genes from young C57BL/6J and aged C57BL/6J datasets were selected, followed by a canonical correlation analysis (CCA) to identify common sources of variation between the two datasets and minimize the batch effect. The first 20 CCA results were chosen for principal component analysis (PCA). Cells were used for 2-dimensional t-Distributed Stochastic Neighbor Embedding (tSNE) (ref van der maaten and hinton 2008) with 0.8 resolution.

 After running this script, a `mouse_eyes_alignment.Rda` file will be generated inside **_data_** folder.
 Do not modify any files in **_data_** folder.
 
 
### 2. Identify_Cell_Types_Manually.R
All clusters are examed against 122(number may change) CD marker genes.
All cell types are predicted by at least two marker genes with the adjusted p-value smaller than 10^-30.

Endothelial cells were identified by Cdh5, Flt1, Kdr, Pecam1, Plvap, Ptprb, and Vwf.<br />
Pericytes were identified by Dcn, Des, Ifitm1, Mylk, Pdgfrb, and Rgs5.<br />
Hematopoietic cells were identified by Laptm5, Ptprc, and Srgn.<br />
Melanocytes were identified by Mlana and Pmel.<br />
Myelinating Schwann cells were identified by Mbp and Mpz.<br />
Retinal pigment epitheliums were identified by Rlbp1 and Rpe65.<br />

Multiple plots and table will be generated, save them if you want. I prefer to keep the original identity of `mouse_eyes_alignment.Rda` intact for further downstream analysis.

### 3. Differential_analysis.R
#### 3.1~3.3 Visualization
`TSNEPlot()`, `SplitDotPlotGG()`,`ggplot()+LabelUR()+LabelLR()` are implemented for visualising differential expressed genes across conditions.

#### 3.4~3.7 Generate csv files with differential expression comparision
`FindAllMarkers.Write()` can split seurat data by conditions(aged vs. young), find All gene Markers differentially expressed between cluster, and generate csv files in **_output_** folder.

Below is a example of `./output/mouse_eyes.aged.csv` file with first 6 rows.


| row.name | p_val | avg_logFC | pct.1 | pct.2 | p_val_adj | cluster  | gene   | 
| ----- | ------ | -------- | ----  | ----- | --------- | ------- | ------|
|  Plvap   |  0 | 3.251091 | 0.950 | 0.239     |     0 | Endothelial Cells  | Plvap
|  Cldn5   |  0 | 2.898564 | 0.858 | 0.071     |     0 | Endothelial Cells  | Cldn5
|  Plpp1   | 0 | 2.794314 | 0.902 | 0.153     |     0 | Endothelial Cells  | Plpp1
|  Egfl7   |  0 | 2.754571 | 0.984 | 0.092    |      0 | Endothelial Cells  | Egfl7
| Igfbp3   |  0 | 2.559527 | 0.725 | 0.053     |     0 | Endothelial Cells | Igfbp3
|    Eng   |  0 | 2.489143 | 0.956 | 0.153     |    0 | Endothelial Cells   | Eng


The results data frame has the following columns :

p_val : p_val (unadjusted) is calculated using likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013) <br />
avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.<br />
pct.1 : The percentage of cells where the gene is detected in the first group.<br />
pct.2 : The percentage of cells where the gene is detected in the second group.<br />
p_val_adj : Adjusted p-value, based on bonferroni correction using all genes in the dataset.<br />
cluster: either cell type name or original identity in `./data/mouse_eyes_alignment.Rda`, will specify in later section.<br />
row.name and gene column are identical.<br />

#### 3.4 Compare differential expression across all major cell types

Perictyes in 129_B6_aged   <——vs——>  all other cells except Perictyes in 129_B6_aged<br />
Endothelial in 129_B6_aged  <——vs——>  all other cells except Endothelial in 129_B6_aged<br />
Myeloid Cells in 129_B6_aged   <——vs——>  all other cells except Myeloid Cells in 129_B6_aged<br />
etc...<br />
Data is stored in `./output/mouse_eyes.aged.csv` and `./output/mouse_eyes.young.csv`.<br />
Cluster names are cell type name.

#### 3.5 Compare differential expression between subcluster within all major cell types

Perictyes subclusters 1 in 129_B6_aged   <——vs——>  all other subclusters in Perictyes except  1 in 129_B6_aged<br />
Endothelial subclusters 1 in 129_B6_aged  <——vs——>  all other subclusters in Endothelial except 1 in 129_B6_aged<br />
Myeloid Cells  subclusters 1 in 129_B6_aged   <——vs——>  all other subclusters in Myeloid Cells except  1 in 129_B6_aged<br />
RPE cells subclusters 1 in 129_B6_aged   <——vs——>  all other subclusters in RPE Cells except  1 in 129_B6_aged<br />
etc...<br />

Respectively, data is stored in :<br />
`./output/Pericytes.aged.csv` and `./output/Pericytes.young.csv`. <br />
`./output/Endothelial.Cells.aged.csv` and `./output/Endothelial.Cells.young.csv`.<br />
`./output/Myeloid.cells.aged.csv` and `./output/Myeloid.cells.young.csv`<br />
`./output/RPE.cells.aged.csv` and `./output/RPE.cells.young.csv`<br />
etc...<br />
Cluster are original identity in `./data/mouse_eyes_alignment.Rda`, from 0 to 16. <br />

Only RPE(Retinal Pigment Epithelium) are further subjected to a second round of unsupervised analysis following the same approach, resulting in 3 subcluster a tSNE analysis with ~0.05 resolution. 

#### 3.6 Compare differential expression in all major cell types across conditions

Perictyes in 129_B6   <——vs——>  Perictyes in 129_B6_aged<br />
Endothelial in 129_B6 <——vs——>  Endothelial in 129_B6_aged<br />
Myeloid Cells in 129_B6   <——vs——>  Myeloid Cells in 129_B6_aged<br />
etc...

#### 3.7 Compare differential expression in all major cell types across conditions

Perictyes subcluster 1 in 129_B6   <——vs——>  Perictyes subcluster 1 in 129_B6_aged<br />
Perictyes subcluster 2 in 129_B6   <——vs——>  Perictyes subcluster 2 in 129_B6_aged<br />
etc...<br />
Endothelial subcluster 1 in 129_B6 <——vs——>  Endothelial subcluster 1 in 129_B6_aged<br />
Endothelial subcluster 2 in 129_B6 <——vs——>  Endothelial subcluster 2 in 129_B6_aged<br />
etc...<br />
and the same for all subclusters in Myeloid cells and RPE cells.
