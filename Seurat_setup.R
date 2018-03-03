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
  colnames(mouse_eyes_raw[[i]]) <- paste0(protocols[i],
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

# lastly, we set the 'protocol' in each dataset for easy identification
# run a canonical correlation analysis to identify common sources
# of variation between the two datasets.
mouse_eyes <- RunCCA(mouse_eyes_Seurat[[1]],mouse_eyes_Seurat[[2]],
                          genes.use = genes.use,
                          num.cc = 30)
save(mouse_eyes, file = "./output/mouse_eyes_alignment.Rda")

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
save(mouse_eyes, file = "./output/mouse_eyes_alignment.Rda")

#====== 1.6 identify phenotype for each cluster  ==========================================
lnames = load(file = "./output/mouse_eyes_alignment.Rda")
lnames

# all marker genes 
marker.genes <- c('Cdh5','Pecam1','Flt1','Vwf','Plvap','Kdr','EMCN',
                  'Car4','ptprb','GPX3','DCN','COL6A1','TIMP3','PDGFRA',
                  'Acta2','Myh11','FGF1','FGF9','SFRP1','PTPRC',
                  'LAPTM5','SRGN','Cma1','Mcpt4','Tpsb2','Cpa3',
                  'MS4A7','LYZ','CD14','CCL2','Emr1','CD68','MARCO','LYZ',
                  'FCGR3A','Itgax','FCER1A','CD3G','Cd4','CD62L',
                  'IL7R','IL2RG','SELL','CREM','Cd8a','CREM','SELL',
                  'Foxp3','Cd19','CD79A','MS4A1','GNLY','Ncr1','KLRD1',
                  'NKG7','Rbfox3','Uty','Ddx3y','Xist','Pdgfrb','Cspg4',
                  'Anpep','Rgs5','Myh11','Mylk','Sost','Des','Vtn',
                  'Ifitm1','Pdgfrb','Vim','Has2','Dcn','Eng','CD44',
                  'Nt5e','Alcam','Kit','Ly6a','Thy1','Itgb1','Vcam1',
                  'Icam2','CD72','CD2','KRT19','Epcam','KRT5','MUC1 ',
                  'SCGB3A2','SCGB1A1','SCGB3A1','SFTPB','FOXJ1 ',
                  'Rpe65','Rlbp1','Msln','Upk3b','Lrrn4','Ccl25',
                  'Cxcl12','Ctsl','Psmb11','Aire','HLA-DMA','Krt5',
                  'Gas1','Plet1','Ly6d','Spink5','Reg3g','Bpifa1',
                  'Pmel','Mlana','MBP','MPZ','SLC36A2','P2RX5','TRIP4',
                  'ASCC1','MYF5','UCP1')
length(marker.genes)
marker.genes <- MouseGenes(mouse_eyes,marker.genes)
#for(i in 1:(length(marker.genes))) print(paste0(i," ",marker.genes[i]))
#Endothelial Cells
FeaturePlot(object = mouse_eyes, 
            features.plot = marker.genes[1:9], min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
#Pericytes and Mesenchymal Cells
pericytes <- marker.genes[c(51,54:56,58:61,64)]
FeaturePlot(object = mouse_eyes, 
            features.plot = pericytes, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# RPE ,Melanocytes and Myelinating Schwann cells
RPE <- marker.genes[c(83:84,92:95)]
FeaturePlot(object = mouse_eyes, 
            features.plot = RPE, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
# Hematopoietic cells
Hematopoietic <- marker.genes[c(19:21,27:28,31,33,39,46)]
FeaturePlot(object = mouse_eyes, 
            features.plot = Hematopoietic, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)


table(mouse_eyes@ident)
old.ident.ids <- 0:16
new.cluster.ids <- c("Retinal Pigment Epithelium",
                     "Pericytes",
                     "Endothelial Cells",
                     "Pericytes",
                     "Pericytes",
                     "Endothelial Cells",
                     "Endothelial Cells",
                     "Pericytes",
                     "Retinal Pigment Epithelium",
                     "Myeloid Cells",
                     "Endothelial Cells",
                     "Myelinating\nSchwann cells",
                     "Myeloid Cells",
                     "Lymphoid Cells",
                     "Pericytes &\n Mesenchymal Cells",
                     "Myelinating\n Schwann cells",
                     "Melanocytes")

mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
# mouse_eyes <- RenameIdentBack(mouse_eyes)
# How many cells are in each cluster
table(mouse_eyes@ident)
TSNEPlot(object = mouse_eyes, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
  ggtitle("Identify conserved cell type markers")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle


# We can also compare proportional shifts in the data. As can be seen in the barplot, 
# the two patients profiled have very different composition
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("Retinal Pigment Epithelium",
                    "Pericytes",
                    "Endothelial Cells",
                    "Myeloid Cells",
                    "Myelinating Schwann cells",
                    "Lymphoid Cells",
                    "Pericytes",
                    "Myelinating Schwann cells",
                    "Melanocytes")
mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)


freq_table <- prop.table(x = table(mouse_eyes@ident, mouse_eyes@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table

#====== 1.7 PCA graphs for marker genes  ==========================================
# Neurons & Pericytes & Mesenchymal Stem Cells & Y-linked & X-linked (only females)
marker.genes <- c("Ihh","Gli1", "Ptch1", "Hhip","Pdgfrb","Cspg4","Anpep","Des","Vtn","Ifitm1")
marker.genes <- MouseGenes(mouse_eyes,marker.genes)

FeaturePlot(object = mouse_eyes, 
            features.plot = marker.genes, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)


#======================================================================
#detect changes in gene expression between 129_B6 and 129_B6_aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to 129_B6_aged or viceversa. 
markers.to.plot <- c("Pmel","Laptm5","Mbp", "Sfrp1","Cd14", "Flt1", "Kdr", "Vwf",
                     "Dcn", "Rgs5","Rpe65")
sdp <- SplitDotPlotGG(mouse_eyes, genes.plot = rev(markers.to.plot),
                      cols.use = c("grey","blue"), x.lab.rot = T, plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")
all.cell <- FetchData(mouse_eyes,"conditions")
cell.129_B6.young <- rownames(all.cell)[all.cell$conditions =="young"]
cell.129_B6.aged <- rownames(all.cell)[all.cell$conditions =="aged"]

mouse_eyes.young <- SubsetData(object = mouse_eyes,
                                cells.use =cell.129_B6.young)
mouse_eyes.aged <- SubsetData(object = mouse_eyes,
                                     cells.use =cell.129_B6.aged)
p1 <- TSNEPlot(mouse_eyes.young, do.label = F, do.return = T,
               no.legend = T, pt.size = 1)+
  ggtitle("Young mouse eyes")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle
p2 <- TSNEPlot(mouse_eyes.aged, do.label = F, do.return = T,
               no.legend = T, pt.size = 1)+
  ggtitle("Aged mouse eyes")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2)
TSNEPlot(mouse_eyes.aged, do.label = F, do.return = T,
         no.legend = F, pt.size = 1)+
  theme(text = element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=10)), #larger legend diagram 
         shape = guide_legend(override.aes = list(size=10))) #larger legend diagram 
  

#========================================================
#2.2 - A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.
table(mouse_eyes.129_B6@ident)
table(mouse_eyes.129_B6_aged@ident)
# We can also compare proportional shifts in the data. As can be seen in the barplot, 
# the two patients profiled have very different composition
freq_table <- prop.table(x = table(mouse_eyes@ident, mouse_eyes@meta.data[, "protocol"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table <- as.data.frame(freq_table)


#2.3 Identify differential expressed genes across conditions
all.cell <- FetchData(mouse_eyes,"protocol")
cell.129_B6 <- rownames(all.cell)[all.cell$protocol !="B6"]
mouse_eyes.129_B6 <- SubsetData(object = mouse_eyes,
                                cells.use =cell.129_B6)
Pericytes <- SubsetData(mouse_eyes.129_B6, ident.use = "Pericytes", subset.raw = T)
Pericytes <- SetAllIdent(Pericytes, id = "protocol")
avg.Pericytes <- log1p(AverageExpression(Pericytes, show.progress = FALSE))
colnames(avg.Pericytes) <- c("B6_129","B6_129_aged")
avg.Pericytes$gene <- rownames(avg.Pericytes)
avg.Pericytes$age <- avg.Pericytes$B6_129_aged - avg.Pericytes$B6_129
head(avg.Pericytes[order(avg.Pericytes$age,decreasing = T),"gene"])

RPE <- SubsetData(mouse_eyes.129_B6, ident.use = "Retinal Pigment Epithelium", subset.raw = T)
RPE <- SetAllIdent(RPE, id = "protocol")
avg.RPE <- log1p(AverageExpression(RPE, show.progress = FALSE))
colnames(avg.RPE) <- c("B6_129","B6_129_aged")
avg.RPE$gene <- rownames(avg.RPE)
avg.RPE$age <- avg.RPE$B6_129_aged - avg.RPE$B6_129
head(avg.RPE[order(avg.RPE$age,decreasing = T),"gene"])


genes.to.label1 = c("Ttr","Trf","Clu","Apod")
genes.to.label2 = c("Ptgds","Rgr")
genes.to.label3 = c("Abhd2","Sult1c1","Ppp1r1b","Trf")
genes.to.label4 = c("Serpine3","Bloc1s1")
p1 <- ggplot(avg.Pericytes, aes(B6_129, B6_129_aged)) + geom_point() + ggtitle("Pericytes")
p1 <- LabelUR(p1, genes = genes.to.label1, avg.Pericytes, 
              adj.u.t = 0.3, adj.u.s = 0.23,text.size = 4)
p1 <- LabelUL(p1, genes = genes.to.label2, avg.Pericytes, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25,text.size = 4)
p2 <- ggplot(avg.RPE, aes(B6_129, B6_129_aged)) + geom_point() + ggtitle("Retinal Pigment Epithelium")
p2 <- LabelUR(p2, genes = genes.to.label3, avg.RPE, 
              adj.u.t = 0.3, adj.u.s = 0.23,text.size = 4)
p2 <- LabelUL(p2, genes = genes.to.label4, avg.RPE, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25,text.size = 4)
plot_grid(p1, p2)

##2018-02-22 ========================================================
#2.4  CVS files with the differential expression between 129_B6 and 129_B6_aged.
#We would need the data for all clusters, as well the subclusters (RPE and hematopoietic cells).
load("./mouse_eyes_alignment.Rda")
table(mouse_eyes@ident)
old.ident.ids <- 0:16
new.cluster.ids <- c("Retinal Pigment Epithelium",
                     "Pericytes",
                     "Endothelial Cells",
                     "Pericytes",
                     "Pericytes",
                     "Endothelial Cells",
                     "Endothelial Cells",
                     "Pericytes",
                     "Retinal Pigment Epithelium",
                     "Myeloid Cells",
                     "Endothelial Cells",
                     "Myelinating Schwann cells",
                     "Myeloid Cells",
                     "Lymphoid Cells",
                     "Pericytes",
                     "Myelinating Schwann cells",
                     "Melanocytes")
mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
table(mouse_eyes@ident)
#mouse_eyes <- RenameIdentBack(mouse_eyes)
table(mouse_eyes@ident)
TSNEPlot(object = mouse_eyes, no.legend = FALSE, do.label = TRUE,
         label.size = 6)

#detect changes in gene expression between 129_B6 and 129_B6_aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to 129_B6_aged or viceversa. 
# by email :Thursday, February 22, 2018 at 4:19 PM

FindBothMarkers <- function(object = mouse_eyes){
  all.cell <- FetchData(object,"protocol")
  cell.young <- rownames(all.cell)[all.cell$protocol =="129_B6"]
  cell.aged <- rownames(all.cell)[all.cell$protocol =="129_B6_aged"]
  
  object.young <- SubsetData(object = object,
                             cells.use =cell.young)
  object.aged <- SubsetData(object = object,
                            cells.use =cell.aged)
  object.young.markers <- FindAllMarkers(object = object.young,
                                         thresh.use = -Inf,
                                         test.use = "bimod",
                                         min.pct = -Inf,
                                         min.diff.pct = -Inf,
                                         min.cells = -Inf)
  object.aged.markers <- FindAllMarkers(object = object.aged, 
                                        thresh.use = -Inf,
                                        test.use = "bimod",
                                        min.pct = -Inf,
                                        min.diff.pct = -Inf,
                                        min.cells = -Inf)
  object.markers <- list(young = object.young.markers,
                         aged = object.aged.markers)
  mapply(write.csv,
         x= object.markers,
         #convert variable (object) name into String
         file=paste(deparse(substitute(mouse_eyes)), 
                    names(object.markers),
                    "csv", sep="."))
  return(object.markers)
}
mouse_eyes.129_B6.markers <- FindBothMarkers(object = mouse_eyes)

#http://satijalab.org/seurat/de_vignette.html#perform-de-analysis-using-alternative-tests
# Compare subclusters within aged and young
# Myeloid Cells===============
load("./mouse_eyes_alignment.Rda")
table(mouse_eyes@ident)
TSNEPlot(object = mouse_eyes, no.legend = FALSE, do.label = TRUE,
         label.size = 6)
Myeloid.cells <- SubsetData(object = mouse_eyes, ident.use = c(10,12))
table(Myeloid.cells@ident)
Myeloid.cells.markers <- FindBothMarkers(Myeloid.cells)

# Perictyes==========
Pericytes <- SubsetData(object = mouse_eyes, ident.use = c(0,3,4,6))
table(Pericytes@ident)
Pericytes.markers <- FindBothMarkers(Pericytes)

# Endothelial Cells==========
Endothelial.Cells <- SubsetData(object = mouse_eyes, ident.use = c(2,5,8,11))
table(Endothelial.Cells@ident)
Endothelial.Cells.markers <- FindBothMarkers(Endothelial.Cells)

# RPE cells=========
RPE.cells <- SubsetData(object = mouse_eyes,
                        ident.use = c(1,7))
RPE.cells <- FindClusters(object = RPE.cells, 
                          reduction.type = "tsne", 
                          dims.use = 1:2, 
                          resolution = 0.05,
                          save.SNN = TRUE)
TSNEPlot(object = RPE.cells, no.legend = TRUE, do.label = TRUE,
         label.size = 7)
RPE.cells.markers <- FindBothMarkers(RPE.cells)

# Compare subcluster between aged vs young
