# scRNAseq-MouseEyes

## Research Background
The choroid is a highly vascularized layer of the eye localized between the sclera and the outermost retinal layer, the retinal pigment epithelium (RPE). Blood supplied by choroidal circulation is the main source of oxygen and nutrients for the RPE and photoreceptors, as well as the main evacuation route for retinal waste. On the other hand, the RPE provides essential support functions for photoreceptor homeostasis and visual function. Thus, it is not surprising that RPE/choroid alterations are present in a vast range of posterior ophthalmic diseases including uveal inflammation and macular degenerative diseases such as Best disease, Stargardt disease and age-related macular degeneration (AMD).

AMD affects 25% people over 80 years old, being a leading cause of blindness in developed countries and approaching epidemic proportions in the United States. However, currently, there is the only treatment for 10% of patients who suffer the advanced, angiogenic form of the disease, and not all respond successfully. The lack of therapeutic options is mainly due to the fact that AMD etiology remains unknown, and probably involves the malfunction of several cell types and multiple intercellular signaling pathways. Thus, it is imperative to carry out systematic studies to characterize in detail all RPE/choroid cell types at the molecular level and to understand how cells communicate among them to maintain choroid and retinal homeostasis. This would constitute the first step to assess whether choroidal intercellular cellular crosstalk is compromised in AMD and to develop new therapeutic strategies based on the restoration of such signaling circuits.

Here, we used young and aged mice to carry out for the first time single cell RNAseq of RPE/choroid tissue, providing its molecular characterization at single cell resolution. In combination with tissue-specific endothelial cell transcription profiling and the use of transgenic mouse models, we discovered some novel mechanisms.

## Data
Chromium single-cell RNA-seq outputs were processed by Cell Ranger analysis pipelines. Data is currently unavailable to the public before publication.

## How to use this repository

#### Software Setup
R version 3.4.3 http://cran.us.r-project.org/bin/macosx/R-3.4.3.pkg <br />
dplyr_0.7.4 <br />
Seurat_2.1.0 https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_2.1.0.tar.gz (other Seurat versions will generate slightly different results )<br />

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

#### 3.4~3.7 Generate CSV files with differential expression comparison
`FindAllMarkersInSameAge()` can split Seurat data by age (aged vs. young), find All gene Markers differentially expressed among clusters within the same age, calculate average UMI, and generate CSV files in **_output_** folder.

`FindAllMarkersbyAge()` can rename cell identity by age (aged vs. young), find All gene Markers
differentially expressed between aged and young strains, calculate average UMI, and generate CSV files in **_output_** folder.

Below is a example of `./output/mouse_eyes.aged.csv` file with first 6 rows.

| row.name | p_val | avg_logFC | pct.1 | pct.2 | p_val_adj |  avg_UMI | cluster | gene |
| ----- | ------ | -------- | ----  | ----- | ------- | ------- | ------| --- |
|  Plvap  |   0  | 3.2510 | 0.950 | 0.239  |    0 | 3.3911 | Endothelial Cells | Plvap
|  Cldn5  |   0  | 2.8985 | 0.858 | 0.071   |    0 | 2.5149 | Endothelial Cells | Cldn5
|  Plpp1  |   0  | 2.7943 | 0.902 | 0.153   |    0 | 2.6818 | Endothelial Cells | Plpp1
|  Egfl7  |   0  | 2.7545 | 0.984 | 0.092   |    0 | 2.8285 | Endothelial Cells | Egfl7
| Igfbp3  |   0  | 2.5595 | 0.725 | 0.053   |    0 | 1.8029 | Endothelial Cells | Igfbp3
|    Eng  |   0  | 2.4891 | 0.956 | 0.153   |    0 | 2.6239 | Endothelial Cells  | Eng

The results data frame has the following columns :

p_val: p_val (unadjusted) is calculated using likelihood-ratio test for single-cell gene expression, (McDavid et al., Bioinformatics, 2013) <br />
avg_logFC: log fold-change of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group.<br />
pct.1: The percentage of cells where the gene is detected in the first group.<br />
pct.2: The percentage of cells where the gene is detected in the second group.<br />
p_val_adj: Adjusted p-value, based on Bonferroni correction using all genes in the dataset.<br />
avg_UMI: average UMI of the cluster.<br />
cluster : either cell types or original clusters in `./data/mouse_eyes_alignment.Rda`. Will be specified in later section.<br />
row.name and gene column are identical.<br />

#### 3.4 Compare differential expression across all major cell types

Perictyes in 129_B6_aged   <——vs——>  all other cells except Perictyes in 129_B6_aged<br />
Endothelial in 129_B6_aged  <——vs——>  all other cells except Endothelial in 129_B6_aged<br />
Myeloid Cells in 129_B6_aged   <——vs——>  all other cells except Myeloid Cells in 129_B6_aged<br />
etc...<br />
Data is stored in `./output/mouse_eyes.aged.csv` and `./output/mouse_eyes.young.csv`.<br />
Cluster names are cell type name.

#### 3.5 Compare differential expression between subcluster within all major cell types, plus visualization for major cell types.

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


Below is a example of `./output/Myeloid.cells.aged.csv` file with first 6 rows.

| row.name | p_val | avg_logFC | pct.1 | pct.2 | p_val_adj |  avg_UMI | cluster | gene |
| ----- | ------ | -------- | ----  | ----- | --------- | ------- | ------| --- |
|  Ctsb  | 2.020953e-46  | -2.108154  | 0.882  | 1.000  | 3.352760e-42  | 1.7019446 | 9 | Ctsb
|  Apoe  | 1.254068e-41  | -2.906343  | 0.836  | 1.000  | 2.080499e-37  | 2.0895885 | 9 | Apoe
|  Ctsd  | 1.997701e-37  | -2.856376  | 0.518  | 0.984  | 3.314186e-33  | 0.7269393 | 9 | Ctsd
|  C1qb  | 2.712081e-36 |  -2.070283  | 0.173  | 0.984  | 4.499343e-32  | 0.5121927 | 9 | C1qb
| Ms4a7  | 3.403560e-36 | -2.685984  | 0.109  | 0.968  | 5.646506e-32  | 0.1998492 | 9 | Ms4a7
| Rplp0  | 4.885699e-36  |  1.100676  | 1.000  | 0.952  | 8.105375e-32  | 3.6815406 | 9 | Rplp0

Cluster indicates the original cluster 9 in `./data/mouse_eyes_alignment.Rda`.
Only RPE (Retinal Pigment Epithelium) is further subjected to a second round of unsupervised analysis following the same approach, resulting in 3 subclusters a tSNE analysis with ~0.05 resolution. 


#### 3.6 Compare differential expression in all major cell types between young and aged mouse

Perictyes in 129_B6   <——vs——>  Perictyes in 129_B6_aged<br />
Endothelial in 129_B6 <——vs——>  Endothelial in 129_B6_aged<br />
Myeloid Cells in 129_B6   <——vs——>  Myeloid Cells in 129_B6_aged<br />
etc...
Data is stored in `./output/mouse_eyes_young_vs_aged.csv`.<br />


#### 3.7 Compare differential expression in all major cell types between young and aged mouse

Perictyes subcluster 1 in 129_B6   <——vs——>  Perictyes subcluster 1 in 129_B6_aged<br />
Perictyes subcluster 2 in 129_B6   <——vs——>  Perictyes subcluster 2 in 129_B6_aged<br />
etc...<br />
Endothelial subcluster 1 in 129_B6 <——vs——>  Endothelial subcluster 1 in 129_B6_aged<br />
Endothelial subcluster 2 in 129_B6 <——vs——>  Endothelial subcluster 2 in 129_B6_aged<br />
etc...<br />
and the same for all subclusters in Myeloid cells and RPE cells.<br />

Respectively, data is stored in :<br />
`./output/Pericytes_young_vs_aged.csv`, <br />
`./output/Endothelium_young_vs_aged.csv`,<br />
`./output/Myeloid.cells_young_vs_aged.csv`, <br />
`./output/RPE_young_vs_aged.csv`.<br />


Below is a example of `./output/Myeloid.cells_young_vs_aged.csv` file with first 6 rows.

| row.name | p_val | avg_logFC | pct.1 | pct.2 | p_val_adj |  avg_UMI | cluster | gene |
| ----- | ------ | -------- | ----  | ----- | --------- | ------- | ------| --- |
| Gm10116 | 1.28e-15 | 0.4632 | 1.000 | 0.682 | 2.13e-11 | 2.367 | 9_young_vs_aged | Gm10116
|  Tpt1 | 1.62e-15 | 0.6115 | 1.000 | 0.991 | 2.68e-11 | 3.488 | 9_young_vs_aged  |  Tpt1
|  Rplp1 | 2.69e-14 | 0.7171 | 1.000 | 1.000 | 4.46e-10 | 3.517 | 9_young_vs_aged  | Rplp1
|  Rps12 | 9.85e-14 | 0.7310 | 1.000 | 0.982 | 1.63e-09 | 2.791 | 9_young_vs_aged  | Rps12
|   Pfn1 | 4.52e-13 | -0.8179 | 0.944 | 0.973 | 7.50e-09 | 2.305 | 9_young_vs_aged  | Pfn1
|   Ttr | 3.45e-12 | -1.3104 | 0.870 | 1.000 | 5.73e-08 | 2.485 | 9_young_vs_aged  | Ttr


Cluster indicates the original cluster 9 in `./data/mouse_eyes_alignment.Rda`, young vs. aged.<br />
Only RPE (Retinal Pigment Epithelium) is further subjected to a second round of unsupervised analysis following the same approach, resulting in 3 subclusters a tSNE analysis with ~0.05 resolution.

##  
sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.4

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_0.7.4         Seurat_2.1.0        Biobase_2.38.0      BiocGenerics_0.24.0
[5] Matrix_1.2-14       cowplot_0.9.2       ggplot2_2.2.1      

loaded via a namespace (and not attached):
  [1] diffusionMap_1.1-0   Rtsne_0.14           VGAM_1.0-5           colorspace_1.3-2    
  [5] ggridges_0.5.0       class_7.3-14         modeltools_0.2-21    mclust_5.4          
  [9] htmlTable_1.11.2     base64enc_0.1-3      proxy_0.4-22         rstudioapi_0.7      
 [13] DRR_0.0.3            flexmix_2.3-14       lubridate_1.7.4      prodlim_2018.04.18  
 [17] mvtnorm_1.0-7        ranger_0.9.0         codetools_0.2-15     splines_3.4.4       
 [21] R.methodsS3_1.7.1    mnormt_1.5-5         doParallel_1.0.11    robustbase_0.93-0   
 [25] knitr_1.20           tclust_1.3-1         RcppRoll_0.2.2       Formula_1.2-2       
 [29] caret_6.0-79         ica_1.0-1            broom_0.4.4          gridBase_0.4-7      
 [33] ddalpha_1.3.3        cluster_2.0.7-1      kernlab_0.9-26       R.oo_1.22.0         
 [37] sfsmisc_1.1-2        compiler_3.4.4       backports_1.1.2      assertthat_0.2.0    
 [41] lazyeval_0.2.1       lars_1.2             acepack_1.4.1        htmltools_0.3.6     
 [45] tools_3.4.4          bindrcpp_0.2.2       igraph_1.1.0         gtable_0.2.0        
 [49] glue_1.2.0           reshape2_1.4.3       Rcpp_0.12.16         NMF_0.21.0          
 [53] trimcluster_0.1-2    gdata_2.18.0         ape_5.1              nlme_3.1-137        
 [57] iterators_1.0.9      fpc_2.1-11           psych_1.8.3.3        timeDate_3043.102   
 [61] gower_0.1.2          stringr_1.3.0        irlba_2.3.2          rngtools_1.2.4      
 [65] gtools_3.5.0         DEoptimR_1.0-8       MASS_7.3-50          scales_0.5.0        
 [69] ipred_0.9-6          RColorBrewer_1.1-2   yaml_2.1.19          pbapply_1.3-4       
 [73] gridExtra_2.3        pkgmaker_0.22        segmented_0.5-3.0    rpart_4.1-13        
 [77] latticeExtra_0.6-28  stringi_1.1.7        foreach_1.4.4        checkmate_1.8.5     
 [81] caTools_1.17.1       ggjoy_0.4.0          lava_1.6.1           geometry_0.3-6      
 [85] dtw_1.18-1           SDMTools_1.1-221     rlang_0.2.0          pkgconfig_2.0.1     
 [89] prabclus_2.2-6       bitops_1.0-6         lattice_0.20-35      ROCR_1.0-7          
 [93] purrr_0.2.4          bindr_0.1.1          recipes_0.1.2        htmlwidgets_1.2     
 [97] tidyselect_0.2.4     CVST_0.2-1           plyr_1.8.4           magrittr_1.5        
[101] R6_2.2.2             gplots_3.0.1         Hmisc_4.1-1          dimRed_0.1.0        
[105] sn_1.5-2             withr_2.1.2          pillar_1.2.2         foreign_0.8-70      
[109] mixtools_1.1.0       survival_2.42-3      scatterplot3d_0.3-41 abind_1.4-5         
[113] nnet_7.3-12          tsne_0.1-3           tibble_1.4.2         KernSmooth_2.23-15  
[117] grid_3.4.4           data.table_1.10.4-3  FNN_1.1              ModelMetrics_1.1.0  
[121] digest_0.6.15        diptest_0.75-7       xtable_1.8-2         numDeriv_2016.8-1
[125] tidyr_0.8.0          R.utils_2.6.0        stats4_3.4.4         munsell_0.4.3    
[129] registry_0.5         magic_1.5-8  