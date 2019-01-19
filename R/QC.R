########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(kableExtra)
source("../R/Seurat_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
if(!dir.exists("./data/")) dir.create("data")
########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======1.1 Load the data files and Set up Seurat object =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/171011_Single_cell_sample list.xlsx")
colnames(df_samples) <- colnames(df_samples) %>% tolower
sample_n = which(df_samples$tests %in% paste0("test",1))
df_samples[sample_n,] %>% kable() %>% kable_styling()
table(df_samples$tests);nrow(df_samples)
samples <- df_samples$sample[sample_n]
sample.id <- df_samples$sample.id[sample_n]
projects <- df_samples$project[sample_n]

# check missing data
(current <- list.files("data")[!grepl(".Rda|RData",list.files("data"))])
(missing_data <- sample.id[!(sample.id %in% current)])

if(length(missing_data)>0){
        # Move files from Download to ./data and rename them
        species <- "mm10"
        for(missing_dat in missing_data){
                old.pth  <- paste("~/Downloads", missing_dat,"outs",
                                  "filtered_gene_bc_matrices",species,
                                  sep = "/")
                list.of.files <- list.files(old.pth)
                new.folder <- paste("./data", missing_dat,"outs",
                                    "filtered_gene_bc_matrices",
                                    species,sep = "/")
                if(!dir.exists(new.folder)) dir.create(new.folder, recursive = T)
                # copy the files to the new folder
                file.copy(paste(old.pth, list.of.files, sep = "/"), new.folder)
                list.files(new.folder)
        }
}

## Load the dataset
Seurat_raw <- list()
Seurat_list <- list()
species <- "mm10"
for(i in 1:length(samples)){
        Seurat_raw[[i]] <- Read10X(data.dir = paste0("data/",sample.id[i],
                                   "/outs/filtered_gene_bc_matrices/",species))
        colnames(Seurat_raw[[i]]) = paste0(samples[i],"_",colnames(Seurat_raw[[i]]))
        rownames(Seurat_raw[[i]]) = gsub(species,"_",rownames(Seurat_raw[[i]]))
        Seurat_list[[i]] <- CreateSeuratObject(Seurat_raw[[i]],
                                                 project = projects[i],
                                                 min.cells = 200,
                                                 min.genes = 3)
}
remove(Seurat_raw);GC()
#======1.1.2 QC before merge =========================
cell.number <- sapply(Seurat_list, function(x) ncol(x@data))
QC_list <- lapply(Seurat_list, function(x) as.matrix(x@data))
median.nUMI <- sapply(QC_list, function(x) median(colSums(x)))
median.nGene <- sapply(QC_list, function(x) median(apply(x,2,function(y) sum(length(y[y>0])))))

min.nUMI <- sapply(QC_list, function(x) min(colSums(x)))
min.nGene <- sapply(QC_list, function(x) min(apply(x,2,function(y) sum(length(y[y>0])))))

QC.list <- cbind(df_samples[sample_n,],cell.number, median.nUMI, median.nGene, 
                 min.nUMI,min.nGene, row.names = samples)
write.csv(QC.list,paste0(path,"QC_list.csv"))
QC.list %>% kable() %>% kable_styling()

remove(QC_list,median.nUMI,median.nGene,min.nUMI,min.nGene,QC.list);GC()

#========1.1.3 merge ===================================
object <- Reduce(function(x, y) MergeSeurat(x, y, do.normalize = F), Seurat_list)
remove(Seurat_list);GC()

# read and select mitochondial genes
(mito.genes <- grep(pattern = "^mt-", x = rownames(object@raw.data), value = TRUE))
percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ])/Matrix::colSums(object@raw.data)

object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")
object@ident = factor(object@ident,levels = samples)

g1 <- lapply(c("nGene", "nUMI", "percent.mito"), function(features){
        VlnPlot(object = object, features.plot = features, nCol = 3, 
                point.size.use = 0.2,size.x.use = 10, group.by = "ident",
                x.lab.rot = T, do.return = T)
        })
save(g1,file= paste0(path,"g1_2_20180414.Rda"))
