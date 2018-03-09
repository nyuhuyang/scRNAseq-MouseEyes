# combine FindAllMarkers and calculate average UMI

FindAllMarkers.UMI <- function (object, genes.use = NULL, logfc.threshold = 0.25, 
          test.use = "wilcox", min.pct = 0.1, min.diff.pct = -Inf, 
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
            FindMarkers(object = object, assay.type = assay.type, 
                        ident.1 = idents.all[i], ident.2 = NULL, genes.use = genes.use, 
                        logfc.threshold = logfc.threshold, test.use = test.use, 
                        min.pct = min.pct, min.diff.pct = min.diff.pct, 
                        print.bar = print.bar, min.cells = min.cells, 
                        latent.vars = latent.vars, max.cells.per.ident = max.cells.per.ident, 
                        ...)
        }, error = function(cond) {
            return(NULL)
        })
        pct.1_UMI <-rowMeans(as.matrix(x = object@data[, WhichCells(object = object,
                                                                       ident = idents.all[[i]])]))
        pct.2_UMI <-rowMeans(as.matrix(x = object@data[, WhichCells(object = object,
                                                                         ident.remove = idents.all[[i]])]))
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


# find and print differentially expressed genes within all major cell types ================
# combine SubsetData, FindAllMarkers, write.csv parallely

FindAllMarkersInSameAge <- function(object,write.csv = TRUE){
        all.cell <- FetchData(object,"conditions")
        cell.young <- rownames(all.cell)[all.cell$conditions =="young"]
        cell.aged <- rownames(all.cell)[all.cell$conditions =="aged"]
        
        object.young <- SubsetData(object = object,
                                 cells.use =cell.young)
        object.aged <- SubsetData(object = object,
                                cells.use =cell.aged)
        object.young.markers <- FindAllMarkers.UMI(object = object.young,
                                             logfc.threshold = -Inf,
                                             return.thresh = 1,
                                             test.use = "bimod",
                                             min.pct = -Inf,
                                             min.diff.pct = -Inf,
                                             min.cells = -Inf)
        object.aged.markers <- FindAllMarkers.UMI(object = object.aged, 
                                            logfc.threshold = -Inf,
                                            return.thresh = 1,
                                            test.use = "bimod",
                                            min.pct = -Inf,
                                            min.diff.pct = -Inf,
                                            min.cells = -Inf)
        object.markers <- list(young = object.young.markers,
                             aged = object.aged.markers)
        if(write.csv){
          mapply(write.csv,
                 x= object.markers,
                 #convert variable (object) name into String
                 file=paste0("./output/",
                             deparse(substitute(object)), 
                             ".",names(object.markers),
                             ".csv"))}
        return(object.markers)
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
        cells <- FetchData(object,"conditions")
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



# modified GenePlot
# GenePlot.1(do.hover = TRUE) will return ggplot 
GenePlot.1 <- function (object, gene1, gene2, cell.ids = NULL, col.use = NULL, 
                        pch.use = 16, pt.size = 1.5, use.imputed = FALSE, use.scaled = FALSE, 
                        use.raw = FALSE, do.hover = FALSE, data.hover = "ident", 
                        do.identify = FALSE, dark.theme = FALSE, do.spline = FALSE, 
                        spline.span = 0.75, ...) 
{
        cell.ids <- SetIfNull(x = cell.ids, default = object@cell.names)
        data.use <- as.data.frame(x = FetchData(object = object, 
                                                vars.all = c(gene1, gene2), cells.use = cell.ids, use.imputed = use.imputed, 
                                                use.scaled = use.scaled, use.raw = use.raw))
        data.plot <- data.use[cell.ids, c(gene1, gene2)]
        names(x = data.plot) <- c("x", "y")
        ident.use <- as.factor(x = object@ident[cell.ids])
        if (length(x = col.use) > 1) {
                col.use <- col.use[as.numeric(x = ident.use)]
        }
        else {
                col.use <- SetIfNull(x = col.use, default = as.numeric(x = ident.use))
        }
        gene.cor <- round(x = cor(x = data.plot$x, y = data.plot$y), 
                          digits = 2)
        if (dark.theme) {
                par(bg = "black")
                col.use <- sapply(X = col.use, FUN = function(color) ifelse(test = all(col2rgb(color) == 
                                                                                               0), yes = "white", no = color))
                axes = FALSE
                col.lab = "white"
        }
        else {
                axes = TRUE
                col.lab = "black"
        }
        
        if (dark.theme) {
                axis(side = 1, at = NULL, labels = TRUE, col.axis = col.lab, 
                     col = col.lab)
                axis(side = 2, at = NULL, labels = TRUE, col.axis = col.lab, 
                     col = col.lab)
        }
        if (do.spline) {
                spline.fit <- smooth.spline(x = data.plot$x, y = data.plot$y, 
                                            df = 4)
                loess.fit <- loess(formula = y ~ x, data = data.plot, 
                                   span = spline.span)
                points(x = data.plot$x, y = loess.fit$fitted, col = "darkblue")
        }
        if (do.identify | do.hover) {
                p <- ggplot2::ggplot(data = data.plot, mapping = aes(x = x, 
                                                                     y = y))
                p <- p + geom_point(mapping = aes(x = x,y=y, color = col.use), size = pt.size, 
                                    shape = pch.use )
                p <- p + labs(title = gene.cor, x = gene1, y = gene2)
                return(p)
        }
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


# MouseGenes
# turn list of gene character to uniform mouse gene list format
MouseGenes <- function(seurat.object,marker.genes){
        if(missing(seurat.object)) 
          stop("A seurat object must be provided first")
        if(class(seurat.object) != "seurat") 
          stop("A seurat object must be provided first")
        if(missing(marker.genes)) 
          stop("A list of marker genes must be provided")
        # marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
        marker.genes <- as.character(marker.genes)
        marker.genes <- unlist(strsplit(marker.genes,","))
        #        marker.genes <- unique(genes)
        marker.genes <- Hmisc::capitalize(tolower(marker.genes))
        marker.genes <- marker.genes[marker.genes %in% seurat.object@raw.data@Dimnames[[1]]]
        print(length(marker.genes))
        return(marker.genes)
}

# rename ident back to 0,1,2,3...
# convert character ident to numeric
# RenameIdent function will generate error if new.cluster.ids overlap with old.ident.ids
# so if want to rename ident smoothly, plyr::mapvalues is still better
RenameIdentBack <- function(object){
        old.ident.ids <- levels(object@ident)
        old.ident.ids <- sort(old.ident.ids)
        new.cluster.ids <- as.numeric(as.factor(old.ident.ids))
        object@ident <- plyr::mapvalues(x = object@ident,
                                        from = old.ident.ids,
                                        to = new.cluster.ids)
        return(object)
}

# FeaturePlot doesn't return ggplot
# SingleFeaturePlot doesn't take seurat object as input
# modified SingleFeaturePlot, take seurat object and return ggplot
SingleFeaturePlot.1 <- function (object = object, feature = feature, pt.size = 1,
                                 dim.1 = 1, dim.2 = 2,pch.use = 16, cols.use  = c("lightgrey","blue"),
                                 cells.use = NULL,dim.codes, min.cutoff = 0, max.cutoff = Inf,
                                 use.imputed = FALSE, reduction.use = "tsne",no.axes = FALSE, no.legend = TRUE, 
                                 dark.theme = FALSE, do.return = FALSE) 
{
        dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                    slot = "key")
        dim.codes <- paste0(dim.code, c(dim.1, dim.2))
        data.plot <- as.data.frame(GetCellEmbeddings(object = object,
                   reduction.type = reduction.use, dims.use = c(dim.1, 
                        dim.2), cells.use = cells.use))
        x1 <- paste0(dim.code, dim.1)
        x2 <- paste0(dim.code, dim.2)
        data.plot$x <- data.plot[, x1]
        data.plot$y <- data.plot[, x2]
        data.plot$pt.size <- pt.size
        names(x = data.plot) <- c("x", "y")
        data.use <- t(x = FetchData(object = object, vars.all = feature, 
                                    cells.use = cells.use, use.imputed = use.imputed))
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
                                            size = pt.size, shape = pch.use) + scale_color_brewer(palette = cols.use)
                }
                else {
                        p <- p + geom_point(mapping = aes(color = col), 
                                            size = pt.size, shape = pch.use) + scale_color_manual(values = cols.use)
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
                                            size = pt.size, shape = pch.use) + scale_color_gradientn(colors = cols.use, 
                                                                                                     guide = guide_colorbar(title = feature))
                }
        }
        if (no.axes) {
                p <- p + labs(title = feature, x = "", y = "") + theme(axis.line = element_blank(), 
                                                                       axis.text.x = element_blank(), axis.text.y = element_blank(), 
                                                                       axis.ticks = element_blank(), axis.title.x = element_blank(), 
                                                                       axis.title.y = element_blank())
        }
        else {
                p <- p + labs(title = feature, x = dim.codes[1], y = dim.codes[2])
        }
        if (no.legend) {
                p <- p + theme(legend.position = "none")
        }
        if (dark.theme) {
                p <- p + DarkTheme()
        }
        return(p)
}



# split by ages and TSNEPlot
TSNEPlotbyAges <- function(object =object){
    
    all.cell <- FetchData(object,"conditions")
    cell.young <- rownames(all.cell)[all.cell$conditions =="young"]
    cell.aged <- rownames(all.cell)[all.cell$conditions =="aged"]
    
    object.young <- SubsetData(object = object,
                               cells.use =cell.young)
    object.aged <- SubsetData(object = object,
                              cells.use =cell.aged)
    
    p1 <- TSNEPlot(object.young, do.label = F, do.return = T,
                   no.legend = T, pt.size = 1)+
        ggtitle(paste0(deparse(substitute(object)), " in Young mouse"))+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
    p2 <- TSNEPlot(object.aged, do.label = F, do.return = T,
                   no.legend = T, pt.size = 1)+
        ggtitle(paste0(deparse(substitute(object)), " in Aged mouse"))+
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5)) #title in middle
    return(plot_grid(p1, p2, nrow =1))
}



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
        cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
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
