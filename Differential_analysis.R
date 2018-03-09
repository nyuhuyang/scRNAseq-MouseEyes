########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#==3.1 detect changes in gene expression between 129_B6 and 129_B6_aged=====
# Rename ident
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
table(mouse_eyes@ident)
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
                     "Pericytes",
                     "Myelinating Schwann cells",
                     "Melanocytes")
mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
table(mouse_eyes@ident)

#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to 129_B6_aged or viceversa. 
all.cell <- FetchData(mouse_eyes,"conditions")
cell.129_B6.young <- rownames(all.cell)[all.cell$conditions =="young"]
cell.129_B6.aged <- rownames(all.cell)[all.cell$conditions =="aged"]

mouse_eyes.young <- SubsetData(object = mouse_eyes,
                               cells.use =cell.129_B6.young)
mouse_eyes.aged <- SubsetData(object = mouse_eyes,
                              cells.use =cell.129_B6.aged)
new.levels <- levels(mouse_eyes.young@ident)
mouse_eyes@ident <- ordered(mouse_eyes@ident,levels = new.levels)
table(mouse_eyes@ident)

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
p3 <- TSNEPlot(mouse_eyes, do.label = F, do.return = T,
               no.legend = T, pt.size = 1)+
    ggtitle("Young and Aged mouse eyes")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle
plot_grid(p1, p2, p3, nrow =2)

TSNEPlot(mouse_eyes.aged, do.label = F, do.return = T,
         no.legend = F, pt.size = 1)+
  theme(text = element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=10)), #larger legend diagram 
         shape = guide_legend(override.aes = list(size=10))) #larger legend diagram 

#==3.2 SplitDotPlotGG======
# The SplitDotPlotGG function can be useful for viewing conserved cell type markers
# across conditions, showing both the expression level and the percentage of cells
# in a cluster expressing any given gene. 
# Here we plot 2-3 strong marker genes for each of our 13 clusters.
markers.to.plot <- c("Pmel","Dcn", "Laptm5","Mbp", "Sfrp1","Cd14", "Flt1", "Kdr", "Vwf",
                     "Rgs5","Rpe65")
sdp <- SplitDotPlotGG(mouse_eyes, genes.plot = rev(markers.to.plot),
                      cols.use = c("grey","blue"), x.lab.rot = T, plot.legend = T,
                      dot.scale = 8, do.return = T, grouping.var = "conditions")

#===3.3 Identify differential expressed genes by ages =========
Pericytes <- SubsetData(mouse_eyes, ident.use = "Pericytes", subset.raw = T)
Pericytes <- SetAllIdent(Pericytes, id = "conditions")
avg.Pericytes <- log1p(AverageExpression(Pericytes, show.progress = FALSE))
colnames(avg.Pericytes) <- c("aged_mouse","young_mouse")
avg.Pericytes$gene <- rownames(avg.Pericytes)
avg.Pericytes$age <- avg.Pericytes$aged_mouse - avg.Pericytes$young_mouse
head(avg.Pericytes[order(avg.Pericytes$age,decreasing = T),"gene"])

RPE <- SubsetData(mouse_eyes, ident.use = "Retinal Pigment Epithelium", subset.raw = T)
RPE <- SetAllIdent(RPE, id = "conditions")
avg.RPE <- log1p(AverageExpression(RPE, show.progress = FALSE))
colnames(avg.RPE) <- c("aged_mouse","young_mouse")
avg.RPE$gene <- rownames(avg.RPE)
avg.RPE$age <- avg.RPE$aged_mouse - avg.RPE$young_mouse
head(avg.RPE[order(avg.RPE$age,decreasing = T),"gene"])


genes.to.label1 = head(avg.Pericytes[order(avg.Pericytes$age,decreasing = T),"gene"],5)
genes.to.label2 = head(avg.Pericytes[order(avg.Pericytes$age,decreasing = F),"gene"],5)
genes.to.label3 = head(avg.RPE[order(avg.RPE$age,decreasing = T),"gene"],5)
genes.to.label4 = head(avg.RPE[order(avg.RPE$age,decreasing = F),"gene"],5)

p1 <- ggplot(avg.Pericytes, aes(aged_mouse, young_mouse)) + geom_point() + ggtitle("Pericytes")
p1 <- LabelUR(p1, genes = genes.to.label1, avg.Pericytes, 
              adj.u.t = 0.3, adj.u.s = 0.23,text.size = 4)
p1 <- LabelUL(p1, genes = genes.to.label2, avg.Pericytes, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25,text.size = 4)
p2 <- ggplot(avg.RPE, aes(aged_mouse, young_mouse)) + geom_point() + ggtitle("Retinal Pigment Epithelium")
p2 <- LabelUR(p2, genes = genes.to.label3, avg.RPE, 
              adj.u.t = 0.15, adj.u.s = 0.1,text.size = 4)
p2 <- LabelUL(p2, genes = genes.to.label4, avg.RPE, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25,text.size = 4)
plot_grid(p1, p2)

#3.4  Compare DE across all major cell types
#We would need the data for all clusters, as well the subclusters (RPE and hematopoietic cells).
#detect changes in gene expression between 129_B6 and 129_B6_aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to 129_B6_aged or viceversa. 
print("3.4 Compare DE across all major cell types")
mouse_eyes.markers <- FindAllMarkersInSameAge(object = mouse_eyes)

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
# Myeloid Cells===============
Myeloid.cells <- SubsetData(object = mouse_eyes,
                            ident.use = old.ident.ids[(new.cluster.ids %in% "Myeloid Cells")])
Hematopoietic.cells <- SubsetData(object = mouse_eyes,
                                  ident.use = c(old.ident.ids[(new.cluster.ids %in% "Myeloid Cells")],
                                                old.ident.ids[(new.cluster.ids %in% "Lymphoid Cells")]))

table(Myeloid.cells@ident)
TSNEPlotbyAges(Myeloid.cells)

table(Hematopoietic.cells@ident)
TSNEPlotbyAges(Hematopoietic.cells)

Myeloid.cells.markers <- FindAllMarkersInSameAge(Myeloid.cells, write.csv = TRUE)

# Perictyes==========
Pericytes <- SubsetData(object = mouse_eyes,
                        ident.use = old.ident.ids[(new.cluster.ids %in% "Pericytes")])
table(Pericytes@ident)
TSNEPlotbyAges(Pericytes)
Pericytes.markers <- FindAllMarkersInSameAge(Pericytes, write.csv = TRUE)

# Endothelial Cells==========
Endothelial.Cells <- SubsetData(object = mouse_eyes,
                                ident.use = old.ident.ids[(new.cluster.ids %in% "Endothelial Cells")])
table(Endothelial.Cells@ident)
TSNEPlotbyAges(Endothelial.Cells)
Endothelial.Cells.markers <- FindAllMarkersInSameAge(Endothelial.Cells, write.csv = TRUE)



# RPE cells=========
RPE.cells <- SubsetData(object = mouse_eyes,
                        ident.use = old.ident.ids[(new.cluster.ids %in% "Retinal Pigment Epithelium")])
RPE.cells <- FindClusters(object = RPE.cells, 
                          reduction.type = "tsne", 
                          dims.use = 1:2, 
                          resolution = 0.05,
                          save.SNN = TRUE)
table(RPE.cells@ident)
TSNEPlotbyAges(RPE.cells)

RPE.cells.markers <- FindAllMarkersInSameAge(RPE.cells, write.csv = TRUE)

#3.6 Compare subcluster between aged vs young===============
print("3.6 Compare subcluster between aged vs young")
# keep the original ident name intact
# Myeloid Cells===============

Myeloid.gde <- FindAllMarkersbyAge(object = Myeloid.cells)
write.csv(x= Myeloid.gde, file="./output/Myeloid.cells_young_vs_aged.csv")

# Perictyes==========
Pericytes.gde <- FindAllMarkersbyAge(object = Pericytes)
write.csv(x= Pericytes.gde, file="./output/Pericytes_young_vs_aged.csv")

# Endothelial Cells==========
table(Endothelial.Cells@ident)
Endothelial.gde <- FindAllMarkersbyAge(Endothelial.Cells)
write.csv(x= Endothelial.gde, file="./output/Endothelium_young_vs_aged.csv")

# RPE cells=========
RPE.gde <- FindAllMarkersbyAge(RPE.cells)
write.csv(x= RPE.gde, file="./output/RPE_young_vs_aged.csv")
