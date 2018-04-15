library(Matrix)
library(Seurat)
mca.matrix <- readRDS(file = "~/Downloads/MCA/MCA_merged_mat.rds")
mca.metadata <- read.csv("~/Downloads/MCA/MCA_All-batch-removed-assignments.csv", 
                         row.names = 1)
mca <- CreateSeuratObject(raw.data = mca.matrix, meta.data = mca.metadata, project = "MouseCellAtlas")
# Only keep annotated cells
mca <- SubsetData(mca, cells.use = rownames(mca@meta.data[!is.na(mca@meta.data$ClusterID), 
                                                          ]), do.clean = TRUE)
# Leaves us with 242k cells
mca
mca <- NormalizeData(object = mca, normalization.method = "LogNormalize", scale.factor = 10000)
mca <- FindVariableGenes(object = mca, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
hv.genes <- head(rownames(mca@hvg.info), 1000)
mito.genes <- grep(pattern = "^mt-", x = rownames(x = mca@data), value = TRUE)
percent.mito <- Matrix::colSums(mca@raw.data[mito.genes, ])/Matrix::colSums(mca@raw.data)
mca <- AddMetaData(object = mca, metadata = percent.mito, col.name = "percent.mito")
mca <- ScaleData(object = mca, genes.use = hv.genes, display.progress = FALSE, 
                 vars.to.regress = "percent.mito", do.par = TRUE, num.cores = 4)
mca <- RunPCA(object = mca, pc.genes = hv.genes, pcs.compute = 100, do.print = TRUE, 
              pcs.print = 1:5, genes.print = 5)
PCElbowPlot(object = mca, num.pc = 100)
PCHeatmap(mca, pc.use = c(1:3, 70:75), cells.use = 500, do.balanced = TRUE)
mca <- FindClusters(object = mca, reduction.type = "pca", dims.use = 1:75, resolution = 3, 
                    save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)
mca <- RunTSNE(object = mca, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE", 
               nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", 
               fast_tsne_path = "/Users/yah2014/src/FIt-SNE/bin/fast_tsne", 
               max_iter = 2000)
library(cowplot)
p1 <- DimPlot(object = mca, reduction.use = "FItSNE", no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = mca, reduction.use = "FItSNE", no.legend = TRUE, group.by = "Tissue", 
              do.return = TRUE, vector.friendly = TRUE, pt.size = 0.1) + ggtitle("Tissue") + 
        theme(plot.title = element_text(hjust = 0.5))
FeaturePlot(mca, c("S100a9", "Sftpc"), reduction.use = "FItSNE", dark.theme = TRUE, 
            pt.size = 0.1, vector.friendly = TRUE)
plot_grid(p1, p2)
FeaturePlot(mca, c("S100a9", "Sftpc"), reduction.use = "FItSNE", dark.theme = TRUE, 
            pt.size = 0.1, vector.friendly = TRUE)
table(mca@ident)

save(mca, file = "./Dropbox/Public/Olivier/R/scRNAseq-MouseAgedEyes/data/mca.Rda")
