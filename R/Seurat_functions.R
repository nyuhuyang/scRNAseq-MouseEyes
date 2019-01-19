#' a supprting function for SingleFeaturePlot.1 and FeatureHeatmap.1
#' Change ggplot color scale to increase contrast gradient
#' #https://github.com/satijalab/seurat/issues/235
#' @param p ggplot object
#' @param alpha.use Define transparency of points
#' @param gradient.use Change fill and colour gradient values
#' @param scaled.expression.threshold Define lower limit of scaled gene expression level
#' @export p ggplot object
ChangeColorScale <- function(p, alpha.use = 0.8,
                             gradient.use = c("yellow", "red"),
                             scaled.expression.threshold = 0) {
        # Order data by scaled gene expresion level
        # Compute maximum value in gene expression
        if (length(p$data$scaled.expression)>0){                # FeatureHeatmap.1
                p$data <- p$data[order(p$data$scaled.expression),]
                max.scaled.exp <- max(p$data$scaled.expression)
        } else if (length(p$data$gene)>0){                   # SingleFeaturePlot.1 
                p$data <- p$data[order(p$data$gene),] 
                max.scaled.exp <- max(p$data$gene) 
        }
        
        # Define lower limit of scaled gene expression level
        if (!is.null(scaled.expression.threshold)) {
                scaled.expression.threshold <- scaled.expression.threshold
        } else if (is.null(scaled.expression.threshold)) {
                if (length(p$data$scaled.expression)>0){
                        scaled.expression.threshold <- min(p$data$scaled.expression)
                } else if (length(p$data$gene)>0) {
                        scaled.expression.threshold <- min(p$data$gene)
                }
        }
        
        # Fill points using the scaled gene expression levels
        p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
        
        # Define transparency of points
        p$layers[[1]]$mapping$alpha <- alpha.use
        
        # Change fill and colour gradient values
        p = p + guides(colour = FALSE)
        p = p + scale_colour_gradientn(colours = gradient.use, guide = F,
                                       limits = c(scaled.expression.threshold,
                                                  max.scaled.exp),
                                       na.value = "grey") +
                scale_fill_gradientn(colours = gradient.use,
                                     name = expression(atop(Scaled, expression)),
                                     limits = c(scaled.expression.threshold,
                                                max.scaled.exp),
                                     na.value = "grey") +
                scale_alpha_continuous(range = alpha.use, guide = F)
        
        # Return plot
        return(p)
}


#' Combine FindAllMarkers and calculate average UMI
#' Modified Seurat::FindAllMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @export to.return data frame
#' @example FindAllMarkers.UMI(mouse_eyes)
FindAllMarkers.UMI <- function (object, genes.use = NULL, logfc.threshold = 0.25, 
          test.use = "bimod", min.pct = 0.1, min.diff.pct = -Inf, 
          print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
          return.thresh = 0.01, do.print = FALSE, random.seed = 1, 
          min.cells = 3, latent.vars = "nUMI", assay.type = "RNA", 
          ...)
{
    data.1 <- GetAssayData(object = object, assay.type = assay.type, 
                           slot = "data")
    genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.1))
    if ((test.use == "roc") && (return.thresh == 0.01)) {
        return.thresh = 0.7
    }
    idents.all <- sort(x = unique(x = object@ident))
    genes.de <- list()
    avg_UMI <- list()
    for (i in 1:length(x = idents.all)) {
        genes.de[[i]] <- tryCatch({
            FindMarkers.1(object = object, assay.type = assay.type, 
                        ident.1 = idents.all[i], ident.2 = NULL, genes.use = genes.use, 
                        logfc.threshold = logfc.threshold, test.use = test.use, 
                        min.pct = min.pct, min.diff.pct = min.diff.pct, 
                        print.bar = print.bar, min.cells = min.cells, 
                        latent.vars = latent.vars, max.cells.per.ident = max.cells.per.ident, 
                        ...)
        }, error = function(cond) {
            return(NULL)
        })
        pct.1_UMI <-rowMeans(expm1(as.matrix(x = object@data[, WhichCells(object = object,
                                                                       ident = idents.all[[i]])])))
        pct.2_UMI <-rowMeans(expm1(as.matrix(x = object@data[, WhichCells(object = object,
                                                                         ident.remove = idents.all[[i]])])))
        avg_UMI[[i]] <-data.frame(pct.1_UMI, pct.2_UMI)
        genes.de[[i]] <- cbind(genes.de[[i]],
                               avg_UMI[[i]][match(rownames(genes.de[[i]]), 
                                                              rownames(avg_UMI[[i]])),])
        if (do.print) {
            print(paste("Calculating cluster", idents.all[i]))
        }
    }
    gde.all <- data.frame()
    for (i in 1:length(x = idents.all)) {
        if (is.null(x = unlist(x = genes.de[i]))) {
            next
        }
        gde <- genes.de[[i]]
        if (nrow(x = gde) > 0) {
            if (test.use == "roc") {
                gde <- subset(x = gde, subset = (myAUC > return.thresh | 
                                                     myAUC < (1 - return.thresh)))
            }
            else {
                gde <- gde[order(gde$p_val, -gde$avg_logFC), 
                           ]
                gde <- subset(x = gde, subset = p_val < return.thresh)
            }
            if (nrow(x = gde) > 0) {
                gde$cluster <- idents.all[i]
                gde$gene <- rownames(x = gde)
            }
            if (nrow(x = gde) > 0) {
                gde.all <- rbind(gde.all, gde)
            }
        }
    }
    if (only.pos) {
        return(subset(x = gde.all, subset = avg_logFC > 0))
    }
    rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
    return(gde.all)
}


# Modified p_val_adj to method = "BH"
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @export gde.all data frame
FindMarkers.1 <- function (object, ident.1, ident.2 = NULL, genes.use = NULL, 
          logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, 
          min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE, 
          max.cells.per.ident = Inf, random.seed = 1, latent.vars = "nUMI", 
          min.cells = 3, pseudocount.use = 1, assay.type = "RNA", 
          ...) 
{
        data.use <- GetAssayData(object = object, assay.type = assay.type, 
                                 slot = "data")
        genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.use))
        methods.noprefiliter <- c("DESeq2", "zingeR")
        if (test.use %in% methods.noprefiliter) {
                genes.use <- rownames(x = data.use)
                min.diff.pct <- -Inf
                logfc.threshold <- 0
        }
        if (length(x = as.vector(x = ident.1) > 1) && any(as.character(x = ident.1) %in% 
                                                          object@cell.names)) {
                cells.1 <- intersect(x = ident.1, y = object@cell.names)
        }
        else {
                cells.1 <- WhichCells(object = object, ident = ident.1)
        }
        if (length(x = as.vector(x = ident.2) > 1) && any(as.character(x = ident.2) %in% 
                                                          object@cell.names)) {
                cells.2 <- intersect(x = ident.2, y = object@cell.names)
        }
        else {
                if (is.null(x = ident.2)) {
                        cells.2 <- WhichCells(object = object, cells.use = setdiff(object@cell.names, 
                                                                                   cells.1))
                }
                else {
                        cells.2 <- WhichCells(object = object, ident = ident.2)
                }
        }
        cells.2 <- setdiff(x = cells.2, y = cells.1)
        if (length(x = cells.1) == 0) {
                print(paste("Cell group 1 is empty - no cells with identity class", 
                            ident.1))
                return(NULL)
        }
        if (length(x = cells.2) == 0) {
                print(paste("Cell group 2 is empty - no cells with identity class", 
                            ident.2))
                return(NULL)
        }
        if (length(cells.1) < min.cells) {
                stop(paste("Cell group 1 has fewer than", as.character(min.cells), 
                           "cells in identity class", ident.1))
        }
        if (length(cells.2) < min.cells) {
                stop(paste("Cell group 2 has fewer than", as.character(min.cells), 
                           " cells in identity class", ident.2))
        }
        thresh.min <- 0
        data.temp1 <- round(x = apply(X = data.use[genes.use, cells.1, 
                                                   drop = F], MARGIN = 1, FUN = function(x) {
                                                           return(sum(x > thresh.min)/length(x = x))
                                                   }), digits = 3)
        data.temp2 <- round(x = apply(X = data.use[genes.use, cells.2, 
                                                   drop = F], MARGIN = 1, FUN = function(x) {
                                                           return(sum(x > thresh.min)/length(x = x))
                                                   }), digits = 3)
        data.alpha <- cbind(data.temp1, data.temp2)
        colnames(x = data.alpha) <- c("pct.1", "pct.2")
        alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
        names(x = alpha.min) <- rownames(x = data.alpha)
        genes.use <- names(x = which(x = alpha.min > min.pct))
        if (length(x = genes.use) == 0) {
                stop("No genes pass min.pct threshold")
        }
        alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, 
                                        FUN = min)
        genes.use <- names(x = which(x = alpha.min > min.pct & alpha.diff > 
                                             min.diff.pct))
        if (length(x = genes.use) == 0) {
                stop("No genes pass min.diff.pct threshold")
        }
        data.1 <- apply(X = data.use[genes.use, cells.1, drop = F], 
                        MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 
                                                                  pseudocount.use))
        data.2 <- apply(X = data.use[genes.use, cells.2, drop = F], 
                        MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) + 
                                                                  pseudocount.use))
        total.diff <- (data.1 - data.2)
        if (!only.pos) 
                genes.diff <- names(x = which(x = abs(x = total.diff) > 
                                                      logfc.threshold))
        if (only.pos) 
                genes.diff <- names(x = which(x = total.diff > logfc.threshold))
        genes.use <- intersect(x = genes.use, y = genes.diff)
        if (length(x = genes.use) == 0) {
                stop("No genes pass logfc.threshold threshold")
        }
        if (max.cells.per.ident < Inf) {
                set.seed(seed = random.seed)
                if (length(cells.1) > max.cells.per.ident) 
                        cells.1 = sample(x = cells.1, size = max.cells.per.ident)
                if (length(cells.2) > max.cells.per.ident) 
                        cells.2 = sample(x = cells.2, size = max.cells.per.ident)
        }
        if (test.use == "bimod") {
                to.return <- DiffExpTest(object = object, assay.type = assay.type, 
                                         cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                         print.bar = print.bar)
        }
        if (test.use == "roc") {
                to.return <- MarkerTest(object = object, assay.type = assay.type, 
                                        cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                        print.bar = print.bar)
        }
        if (test.use == "t") {
                to.return <- DiffTTest(object = object, assay.type = assay.type, 
                                       cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                       print.bar = print.bar)
        }
        if (test.use == "tobit") {
                to.return <- TobitTest(object = object, assay.type = assay.type, 
                                       cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                       print.bar = print.bar)
        }
        if (test.use == "negbinom") {
                to.return <- NegBinomDETest(object = object, assay.type = assay.type, 
                                            cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                            latent.vars = latent.vars, print.bar = print.bar, 
                                            min.cells = min.cells)
        }
        if (test.use == "poisson") {
                to.return <- PoissonDETest(object = object, assay.type = assay.type, 
                                           cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                           latent.vars = latent.vars, print.bar = print.bar)
        }
        if (test.use == "MAST") {
                to.return <- MASTDETest(object = object, assay.type = assay.type, 
                                        cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                        latent.vars = latent.vars, ...)
        }
        if (test.use == "wilcox") {
                to.return <- WilcoxDETest(object = object, assay.type = assay.type, 
                                          cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                          print.bar = print.bar, ...)
        }
        if (test.use == "DESeq2") {
                to.return <- DESeq2DETest(object = object, assay.type = assay.type, 
                                          cells.1 = cells.1, cells.2 = cells.2, genes.use = genes.use, 
                                          ...)
        }
        to.return[, "avg_logFC"] <- total.diff[rownames(x = to.return)]
        to.return <- cbind(to.return, data.alpha[rownames(x = to.return), 
                                                 ])
        to.return$p_val_adj = p.adjust(p = to.return$p_val, method = "BH", 
                                       n = nrow(GetAssayData(object = object, assay.type = assay.type, 
                                                             slot = "data")))
        if (test.use == "roc") {
                to.return <- to.return[order(-to.return$power, -to.return$avg_logFC), 
                                       ]
        }
        else {
                to.return <- to.return[order(to.return$p_val, -to.return$avg_logFC), 
                                       ]
        }
        if (only.pos) {
                to.return <- subset(x = to.return, subset = avg_logFC > 
                                            0)
        }
        return(to.return)
}


#=====Clean memory======================
# extract from WCGNA
GC <- function()
{
        while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1,
                                                             4]) {
        }
}


# find and print differentially expressed genes across conditions ================
# combine SetIdent,indMarkers and avg_UMI
FindAllMarkersbyAge<- function(object, genes.use = NULL, logfc.threshold = -Inf, 
                              test.use = "bimod", min.pct = -Inf, min.diff.pct = -Inf, 
                              print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
                              return.thresh = 1, do.print = FALSE, random.seed = 1, 
                              min.cells = -Inf, latent.vars = "nUMI", assay.type = "RNA", 
                              ...) 
{   
        cells <- FetchData(object, vars.all = "conditions")
        cells$ident <- object@ident
        cells$new.ident <- paste0(cells$ident,"_",cells$conditions)
        object <- SetIdent(object, cells.use = rownames(cells),
                           ident.use = cells$new.ident)
        idents.all <- sort(x = unique(x = cells$ident))
        genes.de <- list()
        avg_UMI <- list()
        for (i in 1:length(x = idents.all)) {
            genes.de[[i]] <- tryCatch({
                FindMarkers(object = object, 
                            ident.1 = paste0(idents.all[i],"_","young"),
                            ident.2 = paste0(idents.all[i],"_","aged"), 
                            logfc.threshold = logfc.threshold, test.use = test.use, 
                            min.pct = min.pct, min.diff.pct = min.diff.pct, 
                            min.cells = min.cells,...)
            }, error = function(cond) {
                return(NULL)
            })
            pct.1_UMI <-rowMeans(as.matrix(x = object@data[, WhichCells(object = object,
                                            ident = paste0(idents.all[i],"_","young"))]))
            pct.2_UMI <-rowMeans(as.matrix(x = object@data[, WhichCells(object = object,
                                            ident = paste0(idents.all[i],"_","aged"))]))
            avg_UMI[[i]] <-data.frame(pct.1_UMI, pct.2_UMI)
            genes.de[[i]] <- cbind(genes.de[[i]],
                                   avg_UMI[[i]][match(rownames(genes.de[[i]]), 
                                                      rownames(avg_UMI[[i]])),])
        }
        gde.all <- data.frame()
        for (i in 1:length(x = idents.all)) {
            if (is.null(x = unlist(x = genes.de[i]))) {
                next
            }
            gde <- genes.de[[i]]
            if (nrow(x = gde) > 0) {
                if (test.use == "roc") {
                    gde <- subset(x = gde, subset = (myAUC > return.thresh | 
                                                         myAUC < (1 - return.thresh)))
                }
                else {
                    gde <- gde[order(gde$p_val, -gde$avg_logFC), 
                               ]
                    gde <- subset(x = gde, subset = p_val < return.thresh)
                }
                if (nrow(x = gde) > 0) {
                    gde$cluster <- paste0(idents.all[i],"_","young_vs_aged")
                    gde$gene <- rownames(x = gde)
                }
                if (nrow(x = gde) > 0) {
                    gde.all <- rbind(gde.all, gde)
                }
            }
        }
        rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
    #    object.name <- 
    #    write.csv(x= gde.all,
        #convert variable (object) name into String
    #              file=paste0("./output/",deparse(substitute(object)),"_young_vs_aged.csv"))
        return(gde.all)
}


#' Extract ColorHexa from Seurat TSNE plot
#' @param object Seurat object
#' @param cells.use only return ColorHexa of selected cells
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @export colors color vector named by cell ID
#' @example gg_colors(mouse_eyes, cells.use = "129_B6_AAACCTGCATGGTCTA")
gg_colors <- function(object = mouse_eyes, cells.use = NULL, no.legend = TRUE, do.label = TRUE,
                      do.return = FALSE, label.size = 6, gg_title="",...){

        g1 <- Seurat::TSNEPlot(object = object, no.legend = no.legend,
                               do.label = do.label,do.return = do.return,
                               label.size = label.size,...)+
                ggtitle(gg_title)+
                theme(text = element_text(size=15),     #larger text including legend title							
                      plot.title = element_text(hjust = 0.5)) #title in middle
        print(g1)
        g <- ggplot2::ggplot_build(g1)
#        print(unique(g$data[[1]]["colour"]))
        colors <- unlist(g$data[[1]]["colour"])
        
        #select color by cells ID
        cells <- Seurat::WhichCells(object)
        names(colors) <- cells
        if(!is.null(cells.use)) colors <- colors[cells.use]
        
        colors <- as.data.frame(table(colors))
        colors <- colors$colors
        print(colors)
        return(as.character(colors))
}


LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
        for (i in genes) {
                x1 <- exp.mat[i, 1]
                y1 <- exp.mat[i, 2]
                plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                                        label = i, size = text.size)
                plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                                                adj.y.s, yend = y1, size = segment.size)
        }
        return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
        return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                          adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
        return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                          adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}


# clean up the gene names for downstream analysis
# turn list of gene character to uniform mouse gene list format
#' @param object Seurat object
#' @param marker.genes gene names, can be one gene or vector. Must be character
#' @param unique TRUE/FALSE, output unique gene name or not
#' @export marker.genes uniform mouse gene that exsit in object@data
#' @example MouseGenes(mouse_eyes,c("Cdh5","Pecam1","Flt1"))
MouseGenes <- function(object, marker.genes, unique =F){
        # marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
        if(missing(object)) 
                stop("A seurat object must be provided first!")
        if(class(object) != "seurat") 
                stop("A seurat object must be provided first!")
        if(missing(marker.genes)) 
                stop("A list of marker genes must be provided!")
        if(object@var.genes[1]==toupper(object@var.genes[1]))
                stop("This is human genome, use HumanGenes() instead!")
        
        marker.genes <- as.character(marker.genes)
        marker.genes <- unlist(strsplit(marker.genes,","))
        marker.genes <- Hmisc::capitalize(tolower(marker.genes))
        marker.genes <- CaseMatch(search = marker.genes, match = rownames(x = object@raw.data))
        if(unique) marker.genes <- unique(marker.genes)
        print(length(marker.genes))
        return(as.character(marker.genes))
}


# FeaturePlot doesn't return ggplot
# SingleFeaturePlot doesn't take seurat object as input
# modified SingleFeaturePlot, take seurat object and return ggplot
SingleFeaturePlot.1 <- function (object = object, feature = feature, pt.size = 1.0,
                                 dim.1 = 1, dim.2 = 2,pch.use = 16, cols.use  = c("lightgrey","blue"),
                                 gradient.use = c("orangered", "red4"),threshold= NULL,text.size=15,
                                 cells.use = NULL,dim.codes, min.cutoff = 0, max.cutoff = Inf,
                                 use.imputed = FALSE, reduction.use = "tsne",no.axes = FALSE, no.legend = TRUE, 
                                 dark.theme = FALSE,title="", do.return = FALSE,x.lim = NULL,
                                 y.lim = NULL, ...) 
{
        dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                    slot = "key")
        dim.codes <- paste0(dim.code, c(dim.1, dim.2))
        data.plot <- as.data.frame(GetCellEmbeddings(object = object,
                                                     reduction.type = reduction.use, 
                                                     dims.use = c(dim.1,dim.2), 
                                                     cells.use = cells.use))
        x1 <- paste0(dim.code, dim.1)
        x2 <- paste0(dim.code, dim.2)
        data.plot$x <- data.plot[, grep(x1,colnames(data.plot),value = T)]
        data.plot$y <- data.plot[, grep(x2,colnames(data.plot),value = T)]
        data.plot$pt.size <- pt.size
        names(x = data.plot) <- c("x", "y")
        data.use <- t(x = FetchData(object = object, vars.all = feature, 
                                    cells.use = cells.use, use.imputed = use.imputed,...))
        data.gene <- na.omit(object = data.frame(data.use[1,])) # Error in data.use[feature, ] : subscript out of bounds
        min.cutoff <- Seurat:::SetQuantile(cutoff = min.cutoff, data = data.gene)
        max.cutoff <- Seurat:::SetQuantile(cutoff = max.cutoff, data = data.gene)
        data.gene <- sapply(X = data.gene, FUN = function(x) {
                return(ifelse(test = x < min.cutoff, yes = min.cutoff, 
                              no = x))
        })
        data.gene <- sapply(X = data.gene, FUN = function(x) {
                return(ifelse(test = x > max.cutoff, yes = max.cutoff, 
                              no = x))
        })
        data.plot$gene <- data.gene
        
        if (length(x = cols.use) == 1) {
                brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
        }
        else {
                brewer.gran <- length(x = cols.use)
        }
        if (all(data.gene == 0)) {
                data.cut <- 0
        }
        else {
                data.cut <- as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.gene), 
                                                             breaks = brewer.gran)))
        }
        data.plot$col <- as.factor(x = data.cut)
        p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
        if (brewer.gran != 2) {
                if (length(x = cols.use) == 1) {
                        p <- p + geom_point(mapping = aes(color = col), 
                                            size = pt.size, shape = pch.use)
                }
                else {
                        p <- p + geom_point(mapping = aes(color = col), 
                                            size = pt.size, shape = pch.use)
                }
        }
        else {
                if (all(data.plot$gene == data.plot$gene[1])) {
                        warning(paste0("All cells have the same value of ", 
                                       feature, "."))
                        p <- p + geom_point(color = cols.use[1], size = pt.size, 
                                            shape = pch.use)
                }
                else {
                        p <- p + geom_point(mapping = aes(color = gene), 
                                            size = pt.size, shape = pch.use)
                }
        }
        if (no.axes) {
                p <- p + labs(title = feature, x = "", y = "") + theme(axis.line = element_blank(), 
                                                                       axis.text.x = element_blank(), 
                                                                       axis.text.y = element_blank(), 
                                                                       axis.ticks = element_blank(), 
                                                                       axis.title.x = element_blank(), 
                                                                       axis.title.y = element_blank())
        }
        else {
                p <- p + labs(x = dim.codes[1], y = dim.codes[2])
        }
        if (no.legend) {
                p <- p + theme(legend.position = "none")
        }
        if (dark.theme) {
                p <- p + DarkTheme()
        }
        p$data <- p$data[order(p$data$gene),]
        p1 <- ChangeColorScale(p, alpha.use = 1,
                               scaled.expression.threshold = threshold,
                               gradient.use = gradient.use)
        if(!is.null(x.lim)) p1 = p1 + xlim(x.lim)
        if(!is.null(y.lim)) p1 = p1 + ylim(y.lim)
        p1 <- p1 +ggtitle(paste0(title,"\n",feature))+
                theme(text = element_text(size=text.size),     #larger text including legend title							
                      axis.text.x = element_text(size=text.size*0.8),
                      axis.text.y = element_text(size=text.size*0.8),
                      plot.title = element_text(hjust = 0.5,size=text.size*1.5)) #title in middle
        return(p1)
}


#' Split seurat cell names by certein criteria
#' A supporting funtion to SplitSeurat
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @export cell.subsets list of subseted cell names by certein conditions, plus levels to split
#' @example SplitCells(mouse_eyes, split.by = "conditions")
SplitCells <- function(object = mouse_eyes, split.by = "conditions"){

    cell.all <- FetchData(object = object, vars.all = split.by)
    if(class(cell.all[,1]) == "numeric"){
            Levels = paste(c("Express no","Express"), split.by)
            cell.subsets[[1]] <- rownames(cell.all)[cell.all[,split.by] == 0]
            cell.subsets[[2]] <- rownames(cell.all)[cell.all[,split.by] > 0]
            cell.subsets[[3]] <- Levels # record conditions in the last return
    }
    if(class(cell.all[,1]) == "factor"){
            Levels <- levels(cell.all[,1])
            cell.subsets <- lapply(Levels, function(x)
                    rownames(cell.all)[cell.all[,split.by] == x])
            cell.subsets[[length(cell.subsets)+1]] <- Levels # record conditions in the last return
    }
    return(cell.subsets)
}


# find and print differentially expressed genes within all major cell types
# combine SubsetData, FindAllMarkers, write.csv
#' @param object Seurat object
#' @param split.by compatible to vars.all in Seurat::FetchData. If split.by is gene name,
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @export gde.all data frame
#' @example FindAllMarkers.UMI(mouse_eyes)
SplitFindAllMarkers <- function(object, split.by = "conditions", write.csv = TRUE,...){
        
        object.subsets <- SplitCells(object = object, split.by = split.by)
        conditions <- object.subsets[[length(object.subsets)]] # levels of conditions
        
        object.markers <- list()
        
        for(i in 1:length(conditions)){ 
                object.markers[[i]] <- FindAllMarkers.UMI(object = object.subsets[[i]],
                                                          logfc.threshold = -Inf,
                                                          return.thresh = 1,
                                                          test.use = "bimod",
                                                          min.pct = -Inf,
                                                          min.diff.pct = -Inf,
                                                          min.cells = -Inf,...)
                if(write.csv) write.csv(object.markers[[i]], 
                                        file=paste0("./output/",
                                                    deparse(substitute(object)), 
                                                    "_",conditions[i],
                                                    ".csv"))
        }
        return(object.markers)
}


#' Split Seurat object by certein criteria
#' A supporting funtion to SplitTSNEPlot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @export object.subsets list of subseted object by certein conditions,
#' plus levels to be splited
#' @example SplitCells(mouse_eyes, split.by = "conditions")
SplitSeurat <- function(object = object, split.by = "conditions"){

        cell.subsets <- SplitCells(object = object, split.by = split.by)
        
        object.subsets <- list()
        for(i in 1:(length(cell.subsets)-1)){
                object.subsets[[i]] <- SubsetData(object, cells.use =cell.subsets[[i]])
        }
        object.subsets[[i+1]] <- cell.subsets[[i+1]] # record conditions in the last return
        return(object.subsets)
}


#' Split Seurat by certein criteria and make tsne plot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @param select.plots output order, default to NULL. If want to change,use c(2,1) for example
#' @param return.data TRUE/FASLE, return splited ojbect or not.
#' @param ... all other parameters are inherited from Seurat::TSNEPlot()
#' @export p ggplot object from TSNEplot
#' @example SplitTSNEPlot(mouse_eyes, split.by = "conditions")
#' @example SplitTSNEPlot(mouse_eyes, split.by = "Rlbp1", select.plots = c(2,1))
SplitTSNEPlot <- function(object = mouse_eyes, split.by = "conditions",
                          select.plots = NULL, return.plots = FALSE,
                          do.label = T, group.by = "ident", no.legend = TRUE,
                          pt.size = 1,label.size = 5,... ){

    object.subsets <- SplitSeurat(object = object, split.by = split.by)
    levels <- object.subsets[[length(object.subsets)]]
    
    p <- list()
    if(is.null(select.plots)) select.plots <- 1:length(levels)
    for(i in 1:length(select.plots)){ 
        p[[i]] <- TSNEPlot(object = object.subsets[[select.plots[i]]],
                           do.label = do.label, group.by = group.by, 
                           do.return = T, no.legend = no.legend,
                           pt.size = pt.size,label.size = label.size,...)+
            ggtitle(levels[select.plots[i]])+
            theme(text = element_text(size=20),     #larger text including legend title							
                  plot.title = element_text(hjust = 0.5)) #title in middle
    }
    p <- p[lapply(p,length)>0] # remove NULL element
    if(return.plots) return(p) else print(do.call(cowplot::plot_grid, p))
}


#' Generate 3D TSNEplot
#' @param object Seurat object after performing RunTSNE(dim.embed = 3)
#' @param ... all other parameters are inherited from Seurat::TSNEPlot()
#' @export p/p3 ggplot object from TSNEplot
TSNEPlot.3D <- function (object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, dim.3 = 3,
          cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE, 
          cols.use = NULL, group.by = "ident", pt.shape = NULL, do.hover = FALSE, 
          data.hover = "ident", do.identify = FALSE, do.label = FALSE, 
          label.size = 4, no.legend = FALSE, no.axes = FALSE, dark.theme = FALSE, 
          plot.order = NULL, plot.title = NULL, ...) 
{
        embeddings.use = GetDimReduction(object = object, reduction.type = reduction.use, 
                                         slot = "cell.embeddings")
        if (length(x = embeddings.use) == 0) {
                stop(paste(reduction.use, "has not been run for this object yet."))
        }
        if (ncol(x = embeddings.use) < 3) {
                stop(paste(reduction.use, "doesn't have the third dimension.
                           Suggest performing RunTSNE(dim.embed = 3)"))
        }
        cells.use <- Seurat:::SetIfNull(x = cells.use, default = colnames(x = object@data))
        dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                    slot = "key")
        dim.codes <- paste0(dim.code, c(dim.1, dim.2))
        data.plot <- as.data.frame(x = embeddings.use)
        cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
        data.plot <- data.plot[cells.use, dim.codes]
        ident.use <- as.factor(x = object@ident[cells.use])
        if (group.by != "ident") {
                ident.use <- as.factor(x = FetchData(object = object, 
                                                     vars.all = group.by)[cells.use, 1])
        }
        data.plot$ident <- ident.use
        data.plot$x <- data.plot[, dim.codes[1]]
        data.plot$y <- data.plot[, dim.codes[2]]
        data.plot$z <- data.plot[, dim.codes[3]]
        data.plot$pt.size <- pt.size
        if (!is.null(plot.order)) {
                if (any(!plot.order %in% data.plot$ident)) {
                        stop("invalid ident in plot.order")
                }
                plot.order <- rev(c(plot.order, setdiff(unique(data.plot$ident), 
                                                        plot.order)))
                data.plot$ident <- factor(data.plot$ident, levels = plot.order)
                data.plot <- data.plot[order(data.plot$ident), ]
        }
        rgl::open3d()
        
        car::scatter3d(x = data.plot$x, y = data.plot$y, z = data.plot$z)
         
                geom_point(mapping = aes(colour = factor(x = ident)), 
                           size = pt.size)
        if (!is.null(x = pt.shape)) {
                shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use, 
                                                                             1]
                if (is.numeric(shape.val)) {
                        shape.val <- cut(x = shape.val, breaks = 5)
                }
                data.plot[, "pt.shape"] <- shape.val
                p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
                        geom_point(mapping = aes(colour = factor(x = ident), 
                                                 shape = factor(x = pt.shape)), size = pt.size)
        }
        if (!is.null(x = cols.use)) {
                p <- p + scale_colour_manual(values = cols.use)
        }
        p2 <- p + xlab(label = dim.codes[[1]]) + ylab(label = dim.codes[[2]]) + 
                zlab(label = dim.codes[[3]]) +
                scale_size(range = c(pt.size, pt.size))
        p3 <- p2 + SetXAxisGG() + SetYAxisGG()  + SetLegendPointsGG(x = 6) + 
                SetLegendTextGG(x = 12) + no.legend.title + theme_bw() + 
                NoGrid()
        p3 <- p3 + theme(legend.title = element_blank())
        if (!is.null(plot.title)) {
                p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
        }
        if (do.label) {
                centers <- data.plot %>% dplyr::group_by(ident) %>% 
                        summarize(x = median(x = x), y = median(x = y))
                p3 <- p3 + geom_point(data = centers, mapping = aes(x = x, 
                                                                    y = y), size = 0, alpha = 0) + geom_text(data = centers, 
                                                                                                             mapping = aes(label = ident), size = label.size)
        }
        if (dark.theme) {
                p <- p + DarkTheme()
                p3 <- p3 + DarkTheme()
        }
        if (no.legend) {
                p3 <- p3 + theme(legend.position = "none")
        }
        if (no.axes) {
                p3 <- p3 + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                                 axis.text.y = element_blank(), axis.ticks = element_blank(), 
                                 axis.title.x = element_blank(), axis.title.y = element_blank(), 
                                 panel.background = element_blank(), panel.border = element_blank(), 
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                 plot.background = element_blank())
        }
        if (do.identify || do.hover) {
                if (do.bare) {
                        plot.use <- p
                }
                else {
                        plot.use <- p3
                }
                if (do.hover) {
                        if (is.null(x = data.hover)) {
                                features.info <- NULL
                        }
                        else {
                                features.info <- FetchData(object = object, 
                                                           vars.all = data.hover)
                        }
                        return(HoverLocator(plot = plot.use, data.plot = data.plot, 
                                            features.info = features.info, dark.theme = dark.theme))
                }
                else if (do.identify) {
                        return(FeatureLocator(plot = plot.use, data.plot = data.plot, 
                                              dark.theme = dark.theme, ...))
                }
        }
        if (do.return) {
                if (do.bare) {
                        return(p)
                }
                else {
                        return(p3)
                }
        }
        if (do.bare) {
                print(p)
        }
        else {
                print(p3)
        }
}