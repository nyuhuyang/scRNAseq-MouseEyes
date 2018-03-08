library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#====== 2.1 identify phenotype for each cluster  ==========================================
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
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

# The SplitDotPlotGG function can be useful for viewing conserved cell type markers
# across conditions, showing both the expression level and the percentage of cells
# in a cluster expressing any given gene. 
# Here we plot 1-3 strong marker genes for each of our 13 clusters.
markers.to.plot <- c("Pmel","Mlana","Mbp","Mpz","Dcn",
                     "Laptm5","Ptprc","Cd14","Cd3g",
                     "Sfrp1","Flt1", "Kdr", "Vwf",
                     "Rpe65","Rlbp1",
                     "Des","Cd14", 
                     "Rgs5")
markers.to.plot <- MouseGenes(mouse_eyes,markers.to.plot)
markers.to.plot <- unique(markers.to.plot)
sdp <- SplitDotPlotGG(mouse_eyes, genes.plot = rev(markers.to.plot),
                      cols.use = c("grey","blue"), x.lab.rot = T, plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")

# Rename ident
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


#====== 2.2 PCA graphs for marker genes  ==========================================
# Neurons & Pericytes & Mesenchymal Stem Cells & Y-linked & X-linked (only females)
marker.genes <- c("Ihh","Gli1", "Ptch1", "Hhip","Pdgfrb","Cspg4","Anpep","Des","Vtn","Ifitm1")
marker.genes <- MouseGenes(mouse_eyes,marker.genes)

FeaturePlot(object = mouse_eyes, 
            features.plot = marker.genes, min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)



# We can also compare proportional shifts in the data. As can be seen in the barplot, 
# the two patients profiled have very different composition
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids <- idents$Var1
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
                     "Pericytes Cells",
                     "Myelinating Schwann cells",
                     "Melanocytes")
mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)

#=====2.3 - A table with the number of cells of each cluster and subcluster, for both B6 and 129_B6 strains.
freq_table <- prop.table(x = table(mouse_eyes@ident, mouse_eyes@meta.data[, "conditions"]), 
                         margin = 2)
barplot(height = freq_table)

freq_table