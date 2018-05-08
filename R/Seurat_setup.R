########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")
########################################################################
#
#  1 Seurat Alignment 
# 
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# Load the mouse.eyes dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
mouse_eyes_raw <- list()
mouse_eyes_Seurat <- list()
protocols <- c("129_B6","B6")
projects <- c("EC-IB-4698","EC-IB-4698")
conditions <- c("129_B6", "B6")
for(i in 1:length(protocols)){
    mouse_eyes_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                     protocols[i],"/outs/filtered_gene_bc_matrices/mm10/"))
    colnames(mouse_eyes_raw[[i]]) <- paste0(conditions[i],
                                            "_",colnames(mouse_eyes_raw[[i]]))
    mouse_eyes_Seurat[[i]] <- CreateSeuratObject(mouse_eyes_raw[[i]],
                                              min.cells = 3,
                                              min.genes = 200,
                                              project = projects[i],
                                              names.delim = "_")
    mouse_eyes_Seurat[[i]]@meta.data$conditions <- conditions[i]
}
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 200, 
                            high.thresholds = Inf)
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, NormalizeData)
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, ScaleData)
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(mouse_eyes_Seurat, function(x) head(rownames(x@hvg.info), 2000))
genes.use <- unique(c(g[[1]],g[[2]]))
for(i in 1:length(conditions)){
    genes.use <- intersect(genes.use, rownames(mouse_eyes_Seurat[[i]]@scale.data))
}
length(genes.use) # 1/10 of total sample size 11212

#======1.2 Perform a canonical correlation analysis (CCA) =========================
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
mouse_eyes <- RunCCA(mouse_eyes_Seurat[[1]],mouse_eyes_Seurat[[2]],
                     genes.use = genes.use,
                     num.cc = 30)
#save(mouse_eyes, file = "./data/mouse_eyes_alignment.Rda")

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = mouse_eyes, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = mouse_eyes, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = mouse_eyes, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = mouse_eyes, reduction.type = "cca", cells.use = 500, dim.use = c(1:3,11:13), 
           do.balanced = TRUE)

DimHeatmap(object = mouse_eyes, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

#======1.3 align seurat objects =========================
mouse_eyes <- CalcVarExpRatio(object = mouse_eyes, reduction.type = "pca",
                              grouping.var = "conditions", dims.use = 1:13)
mouse_eyes <- SubsetData(mouse_eyes, subset.name = "var.ratio.pca",accept.low = 0.5)
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

mouse_eyes <- AlignSubspace(object = mouse_eyes, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:13)
#Now we can run a single integrated analysis on all cells!
mouse_eyes <- RunTSNE(object = mouse_eyes, reduction.use = "cca.aligned", dims.use = 1:13, 
                      do.fast = TRUE)
mouse_eyes <- FindClusters(object = mouse_eyes, reduction.type = "cca.aligned", dims.use = 1:13, 
                           resolution = 0.8, force.recalc = T, save.SNN = TRUE)

p1 <- TSNEPlot(mouse_eyes, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(mouse_eyes, do.label = F, do.return = T, pt.size = 1)
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)
#dev.off()

#======1.4 QC ==================================
mito.genes <- grep(pattern = "^mt-", x = rownames(x = mouse_eyes@data), value = TRUE)
percent.mito <- Matrix::colSums(mouse_eyes@raw.data[mito.genes, ])/Matrix::colSums(mouse_eyes@raw.data)
mouse_eyes <- AddMetaData(object = mouse_eyes, metadata = percent.mito, col.name = "percent.mito")
mouse_eyes <- ScaleData(object = mouse_eyes, genes.use = genes.use, display.progress = FALSE, 
                        vars.to.regress = "percent.mito")
#Now we can run a single integrated analysis on all cells!
VlnPlot(object = mouse_eyes, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = mouse_eyes, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = mouse_eyes, gene1 = "nUMI", gene2 = "nGene")

mouse_eyes <- FilterCells(object = mouse_eyes, subset.names = c("nGene", "percent.mito"), 
                          low.thresholds = c(600, -Inf), high.thresholds = c(5000, 0.10))

#Now, we annotate the clusters as before based on canonical markers.
#png('./output/TSNEPlot.png')
TSNEPlot(object = mouse_eyes,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
        ggtitle("mouse eyes")+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
save(mouse_eyes, file = "./data/mouse_eyes_alignment.Rda")
