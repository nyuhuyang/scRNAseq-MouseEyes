########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kable)
library(kableExtra)
source("./R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#========Fig 1b ======================
lnames = load(file = "./data/mouse_eyes_CCA_20180414.Rda")
lnames
mouse_eyes <- SetAllIdent(mouse_eyes, id = "ClusterNames_0.8")
table(mouse_eyes@ident)
jpeg(paste0(path,"1B_tSNE.jpeg"), units="in", width=10, height=7,res=600)
TSNEPlot(mouse_eyes,do.label =F)
dev.off()

#========Fig 3a EC ======================
jpeg(paste0(path,"3A_Featureplot.jpeg"), units="in", width=10, height=7,res=600)
FeaturePlot(mouse_eyes, features.plot = c("Ihh","Gli1"),min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 0.5)
dev.off()

EC <- SubsetData(mouse_eyes, ident.use = 2:4)
jpeg(paste0(path,"S3A_RidgePlot_EC.jpeg"), units="in", width=10, height=7,res=600)
RidgePlot(object = EC, features.plot = MouseGenes(EC,c("Ihh","Kdr","PECAM1","CDH5")),
          nCol = 2,use.raw = F)
dev.off()


g <- lapply(MouseGenes(EC,c("Ihh","Kdr","PECAM1","CDH5")), function(gene) {
        cell.use <- EC@cell.names[EC@data[gene,] >0]
        sub.EC = SubsetData(EC, cells.use = cell.use)
        RidgePlot(object =sub.EC, features.plot = MouseGenes(EC,gene))
        })
jpeg(paste0(path,"S3A_RidgePlot_EC_positive.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid, g)
dev.off()

split(row.names(EC@meta.data), EC@meta.data$ClusterNames_0.8) %>% 
        lapply(function(cells_use) {
                temp_data <- SubsetData(EC, cells.use = cells_use)
                table(temp_data@data[c("Gli1"),]>0) #%>% prop.table
        })

split(row.names(EC@meta.data), EC@meta.data$ClusterNames_0.8) %>% 
        lapply(function(cells_use) {
                temp_data <- SubsetData(EC, cells.use = cells_use)
                rowMeans(temp_data@data[MouseGenes(EC,c("Ihh","Kdr","PECAM1","CDH5","Gli1")),])
                })

split(row.names(EC@meta.data), EC@meta.data$ClusterNames_0.8) %>% 
        lapply(function(cells_use) {
                temp_data <- SubsetData(EC, cells_use)
                temp_data@data["Gli1",temp_data@data[c("Gli1"),]>0] %>% mean
        })
mouse_eyes_Markers <- FindAllMarkers.UMI(object = mouse_eyes,test.use = "bimod")
write.csv(x= mouse_eyes_Markers, file= paste0(path,"mouse_eyes_Markers.csv"))

#========Fig 3a stromal ======================
stromal <- SubsetData(mouse_eyes, ident.use = 5:8)
jpeg(paste0(path,"S3A_RidgePlot_stromal.jpeg"), units="in", width=10, height=7,res=600)
RidgePlot(object = stromal, features.plot = MouseGenes(stromal,c("Gli1","Ptch1","Pdgfrb","Col1a1")),
          nCol = 2, y.log = F)
dev.off()

g1 <- lapply(MouseGenes(stromal,c("Gli1","Ptch1","Pdgfrb","Col1a1")), function(gene) {
        cell.use <- stromal@cell.names[stromal@data[gene,] >0]
        sub.stromal = SubsetData(stromal, cells.use = cell.use)
        RidgePlot(object =sub.stromal, features.plot = MouseGenes(stromal,gene))
})
jpeg(paste0(path,"S3A_RidgePlot_stromal_positive.jpeg"), units="in", width=10, height=7,res=600)
do.call(plot_grid, g1)
dev.off()

split(row.names(stromal@meta.data), stromal@meta.data$ClusterNames_0.8) %>% 
        lapply(function(cells_use) {
                temp_data <- SubsetData(stromal, cells.use = cells_use)
                table(temp_data@data[c("Gli1"),]>0) %>% prop.table
        })  

split(row.names(stromal@meta.data), stromal@meta.data$ClusterNames_0.8) %>% 
        lapply(function(cells_use) {
                temp_data <- SubsetData(stromal, cells.use = cells_use)
                rowMeans(temp_data@data[MouseGenes(stromal,c("Gli1","Ptch1","Pdgfrb","Col1a1")),])
        }) 

split(row.names(stromal@meta.data), stromal@meta.data$ClusterNames_0.8) %>% 
        lapply(function(cells_use) {
                temp_data <- SubsetData(stromal, cells_use)
                temp_data@data["Cdh5",temp_data@data[c("Cdh5"),]>0] %>% mean
        })


#========Fig 6a ======================
mouse_eyes_Split <- SplitSeurat(object = mouse_eyes, split.by = "conditions")
mouse_eyes_129_B6 <- mouse_eyes_Split[[1]]
mouse_eyes_129_B6.gde <- FindAllMarkers.UMI(object = mouse_eyes_129_B6)
write.csv(x= mouse_eyes_129_B6.gde, file="./output/129_B6.csv")

#  SubsetData and further split Hematopoietic cells ===============
Hema <- SubsetData(object = mouse_eyes_129_B6,
                   ident.use = "11")
Hema <- FindVariableGenes(object = Hema, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
Hema_hv.genes <- head(rownames(Hema@hvg.info), 1000)
set.seed(42)
Hema <- RunPCA(object = Hema, pc.genes = Hema_hv.genes, pcs.compute = 20, do.print = F, 
              pcs.print = 1:5, genes.print = 5)
#PCElbowPlot(object = Hema)
Hema <- FindClusters(object = Hema, reduction.type = "pca", dims.use = 1:5,
                     force.recalc = T,resolution = 1.3, save.SNN = TRUE)
Hema <- RunTSNE(object = Hema, reduction.use = "pca", dims.use = 1:5, 
               do.fast = TRUE)

jpeg(paste0(path,"1B_tSNE.jpeg"), units="in", width=10, height=7,res=600)
TSNEPlot(object = Hema, no.legend = TRUE, do.label = TRUE, pt.size = 2,
         label.size = 7)
dev.off()
# generate csv file

Hema.gde <- FindAllMarkers.UMI(object = Hema,test.use = "bimod")
write.csv(x= Hema.gde, file="./output/129_B6_Hematopoietic.csv")


#=======================
# a threshold in the plot, so that cells that express Ihh or Gli1 above that threshold are labeled in purple.

mouse_eyes <- SetAllIdent(mouse_eyes,id = "ClusterNames_0.8")
TSNEPlot(mouse_eyes, do.label = T)
RidgePlot(mouse_eyes,features.plot = "Ihh")
for(i in c(0,0.1,1,2)){
        print(i)
        jpeg(paste0(path,"Ihh_tSNE_",i,".jpeg"), units="in", width=10, height=7,res=600)
        g <- SingleFeaturePlot.1(mouse_eyes,feature = "Ihh",threshold = i,gradient.use = c("lightgrey", "blue"),
                                 title = paste0("threshold = ",i))
        print(g)
        dev.off()
}

for(i in c(0.5,1,2)){
        print(i)
        jpeg(paste0(path,"Gli1_tSNE_",i,".jpeg"), units="in", width=10, height=7,res=600)
        g <- SingleFeaturePlot.1(mouse_eyes,feature = "Gli1",threshold = i,gradient.use = c("lightgrey", "blue"),
                                 title = paste0("threshold = ",i))
        print(g)
        dev.off()
}
