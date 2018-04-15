library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames

Featureplot <- function(x){
    p <- FeaturePlot(object = mouse_eyes, 
                     reduction.use = "tsne",
                     features.plot = x, min.cutoff = NA, 
                     cols.use = c("lightgrey","blue"), pt.size = 0.5)
    return(p)
}
Adipocytes <- MouseGenes(mouse_eyes,c("SLC36A2","P2RX5","MYF5","UCP1","TRIP4","ASCC1"))
Endothelium <- MouseGenes(mouse_eyes,c("Cdh5","Pecam1","Flt1","Plvap","Kdr","ptprb",
                                        "Vwf","EMCN","Car4"))
Epithelium <- MouseGenes(mouse_eyes,c("KRT19","Epcam","KRT5",
                                       "MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1","Rpe65",
                                       "Rlbp1","Msln","Upk3b","Lrrn4"))
RPE <- MouseGenes(mouse_eyes,c("Rpe65","Rlbp1"))
Fibroblast <- MouseGenes(mouse_eyes,c("FGF1","FGF9","SFRP1"))
Hematopoietic <- MouseGenes(mouse_eyes,c("PTPRC","LAPTM5","SRGN"))
Myeloid <-  MouseGenes(mouse_eyes,c("PPBP","GNG11","HBA2","HBB","Cma1","Mcpt4","Tpsb2",
                                     "Cpa3","LYZ","S100A9","CD14","CCL2","FCGR3A","MS4A7","VMO1"))
Lymphoid <- MouseGenes(mouse_eyes,c("CD3G","CD3D","CD2","Cd19","CD79A","MS4A1",
                                     "GNLY","Ncr1","CCL5","KLRD1","NKG7"))
Melanocytes <- MouseGenes(mouse_eyes,c("Pmel","Mlana"))
Mesenchymal <- MouseGenes(mouse_eyes,c("Pdgfrb","Vim","Has2","Dcn"))
Myelinating_Schwann_cells <- MouseGenes(mouse_eyes,c("MBP","MPZ"))
Pericytes <- MouseGenes(mouse_eyes,c("Pdgfrb","Cspg4","Anpep","Rgs5",
                                      "Myh11","Mylk","Des","Vtn","Ifitm1"))
Smooth_muscle_cells <- MouseGenes(mouse_eyes,c("Acta2","Myh11"))
Stem_cell <- MouseGenes(mouse_eyes,c("POU5F1","FUT4","CD34","PROM1","ABCG2","Runx1","ATXN1",
                                      "Nes","NCAM","NGFR"))
Stromal_fibroblasts <- MouseGenes(mouse_eyes,c("DCN","COL6A1","TIMP3","PDGFRA"))
Neurons <- MouseGenes(mouse_eyes,c("Ihh","Gli1", "Ptch1", "Hhip"))
# Featureplot
Featureplot(Adipocytes) # Adipocytes
Featureplot(Endothelium) # Endothelial Cells
Featureplot(Epithelium) # Epithelium
Featureplot(c(RPE,Melanocytes,Myelinating_Schwann_cells)) # RPE, Melanocytes, Myelinating Schwann cells
Featureplot(Fibroblast) # Fibroblasts
Featureplot(c(Hematopoietic,Myeloid[7:9],Lymphoid[1:3])) # Hematopoietic cells
Featureplot(Myeloid) # Myeloid cells
Featureplot(Lymphoid) # Lymphoid cells
Featureplot(Mesenchymal) # Mesenchymal cells
Featureplot(Pericytes) # Pericytes
Featureplot(Smooth_muscle_cells)
Featureplot(Stem_cell)
Featureplot(Stromal_fibroblasts)
Featureplot(Neurons)

markers.to.plot <- c(Hematopoietic[1:2], Lymphoid[1:2],Myeloid[c(7,9)],
                     Endothelium[c(1:3,5,7)], Myelinating_Schwann_cells,Melanocytes,
                     RPE,Mesenchymal[c(1,4)],Pericytes[c(4,6:7)],
                     Smooth_muscle_cells)
markers.to.plot <- unique(markers.to.plot)
DotPlot(mouse_eyes, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

# Rename ident
table(mouse_eyes@ident)
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("RPE 0",
                     "Mesenchymal cells 1",
                     "Endothelial cells 2",
                     "Endothelial cells 3",
                     "Pericytes 4",
                     "Smooth muscle cells 5",
                     "Mesenchymal cells 6",
                     "Monocytes 7",
                     "RPE 8",
                     "Schwann cells 9",
                     "Endothelial cells 10",
                     "Monocytes 11",
                     "T cells 12")

mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
DotPlot(mouse_eyes, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T,
        dot.scale = 8, do.return = T)
# mouse_eyes <- RenameIdentBack(mouse_eyes)

#====== 2.2 Plots ==========================================
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
table(mouse_eyes@ident)
new.cluster.ids <- c("Retinal Pigment Epithelium",
                     "Mesenchymal cells",
                     "Endothelial cells",
                     "Endothelial cells",
                     "Pericytes",
                     "Smooth\n muscle cells",
                     "Mesenchymal cells",
                     "Monocytes",
                     "Retinal Pigment Epithelium",
                     "Schwann cells",
                     "Endothelial cells",
                     "Monocytes",
                     "T cells")

mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
DotPlot(mouse_eyes, genes.plot = rev(markers.to.plot),
        cols.use = c("blue","red"), x.lab.rot = T, plot.legend = F,
        dot.scale = 8, do.return = T)

TSNEPlot(object = mouse_eyes, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 6)+
  ggtitle("TSNE plot of major cell types")+
  theme(text = element_text(size=20),     #larger text including legend title							
        plot.title = element_text(hjust = 0.5)) #title in middle
freq_table <- prop.table(x = table(mouse_eyes@ident, mouse_eyes@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table

#====== 2.3 Compare cell type changes across conditions  ==========================================
# the two patients profiled have very different composition
# Compare clusters for each dataset
SplitTSNEPlot(mouse_eyes, "conditions")
