########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters (RPE and hematopoietic cells).
#detect changes in gene expression between 129_B6 and B6, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to B6 or viceversa. 

# 3.1.1 load data
# Rename ident
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
table(mouse_eyes@ident)
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("0) Pericytes",
                     "1) Pericytes",
                     "2) Endothelial cells",
                     "3) Smooth muscle cells",
                     "4) Retinal pigment epithelium",
                     "5) Endothelial cells",
                     "6) Pericytes",
                     "7) Pericytes",
                     "8) Hematopoietic cells",
                     "9) Myelinating schwann cells",
                     "10) Endothelial cells",
                     "11) Myelinating schwann cells",
                     "12) Melanocytes")
mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
print("3.1 Compare DE across all major cell types")

# 3.1.2 FindAllMarkers.UMI
mouse_eyes_Split <- SplitCells(object = mouse_eyes, split.by = "conditions")
mouse_eyes_129_B6 <- mouse_eyes_Split[[1]]
mouse_eyes_129_B6.gde <- FindAllMarkers.UMI(object = mouse_eyes_129_B6)
write.csv(x= mouse_eyes_129_B6.gde, file="./output/129_B6.csv")

# 3.2 Compare differential expression between subcluster within all major cell types
# plus visualize all major cell types.
#http://satijalab.org/seurat/de_vignette.html#perform-de-analysis-using-alternative-tests
# Compare subclusters within aged and young
# keep the original ident name intact
print("3.2 Compare DE between subcluster within all major cell types, and visualize all major cell types")
# 3.2.1 SubsetData and further split RPE ===============
RPE <- SubsetData(object = mouse_eyes_129_B6,random.seed = 1,
                  ident.use = "4) Retinal pigment epithelium")
#deselect.cells <- TSNEPlot(object = RPE, do.identify = T)
set.seed(1)
select.cells <- WhichCells(object = RPE)
deselect.cells <- which(select.cells == "129_B6_CGCTATCTCAGGTTCA")
select.cells <- select.cells[-deselect.cells]
RPE <- SubsetData(object = RPE, cells.use =select.cells)
RPE <- FindVariableGenes(object = RPE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
RPE_hv.genes <- head(rownames(RPE@hvg.info), 1000)
RPE <- RunPCA(object = RPE, pc.genes = RPE_hv.genes, pcs.compute = 20, 
              do.print = F)
#PCElbowPlot(object = RPE)
RPE <- FindClusters(object = RPE, reduction.type = "pca", dims.use = 1:5,
                    force.recalc = T, resolution = 0.3, save.SNN = TRUE)
RPE <- RunTSNE(object = RPE, reduction.use = "pca", dims.use = 1:5, 
                      do.fast = TRUE)
TSNEPlot(object = RPE, no.legend = TRUE, do.label = TRUE, pt.size = 2,
         label.size = 7)
# generate csv file
RPE.gde <- FindAllMarkers.UMI(object = RPE,test.use = "bimod")
write.csv(x= RPE.gde, file="./output/129_B6_RPE.csv")

# 3.2.2 SubsetData and further split Hematopoietic cells ===============
Hema <- SubsetData(object = mouse_eyes_129_B6,
                   ident.use = "8) Hematopoietic cells")
Hema <- FindVariableGenes(object = Hema, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
Hema_hv.genes <- head(rownames(Hema@hvg.info), 1000)
set.seed(1)
Hema <- RunPCA(object = Hema, pc.genes = Hema_hv.genes, pcs.compute = 20, do.print = F, 
              pcs.print = 1:5, genes.print = 5)
#PCElbowPlot(object = Hema)
Hema <- FindClusters(object = Hema, reduction.type = "pca", dims.use = 1:5,
                     force.recalc = T,resolution = 1.3, save.SNN = TRUE)
Hema <- RunTSNE(object = Hema, reduction.use = "pca", dims.use = 1:5, 
               do.fast = TRUE)
TSNEPlot(object = Hema, no.legend = TRUE, do.label = TRUE, pt.size = 2,
         label.size = 7)
# generate csv file
Hema.gde <- FindAllMarkers.UMI(object = Hema,test.use = "bimod")
write.csv(x= Hema.gde, file="./output/129_B6_Hematopoietic.csv")