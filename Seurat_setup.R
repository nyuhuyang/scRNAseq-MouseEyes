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
protocols <- c("129_B6","129_B6_aged")
projects <- c("EC-IB-4698","EC-IB-4867")
conditions <- c("young", "aged")
for(i in 1:length(protocols)){
    mouse_eyes_raw[[i]] <- Read10X(data.dir = paste0("./data/",
                                                     protocols[i],"/outs/filtered_gene_bc_matrices/mm10/"))
    colnames(mouse_eyes_raw[[i]]) <- paste0(conditions[i],
                                            "_",colnames(mouse_eyes_raw[[i]]))
}
mouse_eyes_Seurat <- lapply(mouse_eyes_raw, CreateSeuratObject,
                            min.cells = 3,
                            min.genes = 200,
                            project = projects)
for(i in 1:length(protocols)) mouse_eyes_Seurat[[i]]@meta.data$conditions <- conditions[i]
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, FilterCells, 
                            subset.names = "nGene", 
                            low.thresholds = 500, 
                            high.thresholds = Inf)
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, NormalizeData)
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, ScaleData)
mouse_eyes_Seurat <- lapply(mouse_eyes_Seurat, FindVariableGenes, do.plot = FALSE)

# we will take the union of the top 1k variable genes in each dataset for
# alignment note that we use 1k genes in the manuscript examples, you can
# try this here with negligible changes to the overall results
g <- lapply(mouse_eyes_Seurat, function(x) head(rownames(x@hvg.info), 1000))
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
save(mouse_eyes, file = "./data/mouse_eyes_alignment.Rda")

# CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = mouse_eyes, reduction.use = "cca", group.by = "conditions", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = mouse_eyes, features.plot = "CC1", group.by = "conditions", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = mouse_eyes, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

p3 <- MetageneBicorPlot(mouse_eyes, grouping.var = "conditions", dims.eval = 1:30, 
                        display.progress = FALSE)
p3 + geom_smooth(method = 'loess')
DimHeatmap(object = mouse_eyes, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
           do.balanced = TRUE)

DimHeatmap(object = mouse_eyes, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)

PrintDim(object = mouse_eyes, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)


#======1.3 QC (skip)==================================


#======1.4 align seurat objects =========================
#Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned

mouse_eyes <- AlignSubspace(object = mouse_eyes, reduction.type = "cca", grouping.var = "conditions", 
                            dims.align = 1:20)
#Now we can run a single integrated analysis on all cells!
mouse_eyes <- RunTSNE(object = mouse_eyes, reduction.use = "cca.aligned", dims.use = 1:20, 
                      resolution = 0.8, do.fast = TRUE)
mouse_eyes <- FindClusters(object = mouse_eyes, reduction.type = "cca.aligned", dims.use = 1:20, 
                           force.recalc = T, save.SNN = TRUE)
p1 <- TSNEPlot(mouse_eyes, do.return = T, pt.size = 1, group.by = "conditions")
p2 <- TSNEPlot(mouse_eyes, do.label = F, do.return = T, pt.size = 1)
#png('./output/TSNESplot_alignment.png')
plot_grid(p1, p2)
#dev.off()

#Now, we annotate the clusters as before based on canonical markers.
#png('./output/TSNEPlot.png')
TSNEPlot(object = mouse_eyes,do.label = TRUE, group.by = "ident", 
         do.return = TRUE, no.legend = TRUE,
         pt.size = 1,label.size = 8 )+
    ggtitle("TSNE plot of all clusters")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle
#dev.off()
save(mouse_eyes, file = "./data/mouse_eyes_alignment.Rda")