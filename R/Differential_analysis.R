########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#==3.1 detect changes in gene expression between 129_B6 and 129_B6_aged sub celltype=====
# Rename ident
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
table(mouse_eyes@ident)
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Retinal Pigment Epithelium",
                     "Mesenchymal cells",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Pericytes",
                     "Smooth muscle cells",
                     "Mesenchymal cells",
                     "Monocytes",
                     "Retinal Pigment Epithelium",
                     "Schwann cells",
                     "Endothelial cells",
                     "Monocytes",
                     "T cells")

#3.4  Compare DE across all major cell types
#We would need the data for all clusters, as well the subclusters (RPE and hematopoietic cells).
#detect changes in gene expression between 129_B6 and 129_B6_aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to 129_B6_aged or viceversa. 
print("3.4 Compare DE across all major cell types")
mouse_eyes.markers <- SplitFindAllMarkers(object = mouse_eyes)

mouse_eyes.gde <- FindAllMarkersbyAge(object = mouse_eyes)
write.csv(x= mouse_eyes.gde, file="./output/mouse_eyes_young_vs_aged.csv")

# 3.5 Compare differential expression between subcluster within all major cell types
# plus visualize all major cell types.
#http://satijalab.org/seurat/de_vignette.html#perform-de-analysis-using-alternative-tests
# Compare subclusters within aged and young
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
# keep the original ident name intact
print("3.5 Compare DE between subcluster within all major cell types, and visualize all major cell types")
# SubsetData===============
T.cells <- SubsetData(object = mouse_eyes,
                            ident.use = old.ident.ids[(new.cluster.ids %in% "T cells")])
Monocytes <- SubsetData(object = mouse_eyes,
                                  ident.use = c(old.ident.ids[(new.cluster.ids %in% "Monocytes")]))
EC <- SubsetData(object = mouse_eyes,
                 ident.use = old.ident.ids[(new.cluster.ids %in% "Endothelial cells")])
RPE <- SubsetData(object = mouse_eyes,
                 ident.use = old.ident.ids[(new.cluster.ids %in% "Retinal Pigment Epithelium")])

Pericytes <- SubsetData(object = mouse_eyes,
                  ident.use = old.ident.ids[(new.cluster.ids %in% "Mesenchymal cells") | 
                                            (new.cluster.ids %in% "Pericytes") |
                                            (new.cluster.ids %in% "Smooth muscle cells")])
# split TSNE plot=============
SplitTSNEPlot(T.cells)
SplitTSNEPlot(Monocytes)
SplitTSNEPlot(EC)
SplitTSNEPlot(RPE)
SplitTSNEPlot(Pericytes)


# split and find all markers
table(T.cells@ident)
T.cells.markers <- SplitFindAllMarkers(T.cells, write.csv = TRUE)

table(Monocytes@ident)
Monocytes.markers <- SplitFindAllMarkers(Monocytes, write.csv = TRUE)

table(EC@ident)
EC.markers <- SplitFindAllMarkers(EC, write.csv = TRUE)

table(RPE@ident)
RPE.markers <- SplitFindAllMarkers(RPE, write.csv = TRUE)

table(Pericytes@ident)
Pericytes.markers <- SplitFindAllMarkers(Pericytes, write.csv = TRUE)

#3.6 Compare subcluster between aged vs young===============
print("3.6 Compare subcluster between aged vs young")
# keep the original ident name intact
# Myeloid Cells===============

Myeloid.gde <- FindAllMarkersbyAge(object = Monocytes)
write.csv(x= Myeloid.gde, file="./output/Monocytes_young_vs_aged.csv")

# Perictyes==========
Pericytes.gde <- FindAllMarkersbyAge(object = Pericytes)
write.csv(x= Pericytes.gde, file="./output/Pericytes_young_vs_aged.csv")

# Endothelial Cells==========
table(EC@ident)
Endothelial.gde <- FindAllMarkersbyAge(EC)
write.csv(x= Endothelial.gde, file="./output/Endothelium_young_vs_aged.csv")

# RPE cells=========
RPE.gde <- FindAllMarkersbyAge(RPE)
write.csv(x= RPE.gde, file="./output/RPE_young_vs_aged.csv")
