#devtools::load_all("/Users/ekol-yal/Applications/CellAnnotatoRSeurat")

#usethis::use_package(Matrix, lib.loc="/System/Volumes/Data/Users/ekol-yal/opt/anaconda3/pkgs/r-matrix-1.5_1-r42hce01bf1_0/lib/R/library")
usethis::use_package('Seurat')
usethis::use_package('dplyr')
usethis::use_package('ggplot2')

#' Relabel a classification with newly defined labels
#' @param obj a seurat object
#' @param clr2change classification to be relabeled
#' @param new.cluster.ids new label, make sure the order in new.cluster.ids should correspond to that in levels(obj@meta.data[, clr2change])
#' @export changeLabel
changeLabel <- function(obj, clr2change, new.cluster.ids) {
  Idents(obj) <- obj@meta.data[, clr2change]
  names(new.cluster.ids) <- levels(obj@meta.data[, clr2change])
  obj <- Seurat::RenameIdents(obj, new.cluster.ids)
  obj@meta.data[, clr2change] <- Seurat::Idents(obj)
  return(obj)
}


#' Annotate clusters based on markers
#' @param obj a seurat object
#' @param marker.list a list of markers for different cell types
#' @param classification classification to be annotated
#' @export annOnMarker
annOnMarker <- function(obj, marker.list, classification="SCT_snn_res.0.3") {
  idf.wei <- log(1+ncol(obj)/(rowSums(as.matrix(obj@assays$SCT@data>0)) + 1))
  tf.idf <- Matrix(t(as.matrix(obj@assays$SCT@data)), sparse=TRUE)
  #tf.idf <- t(as.matrix(SCT@data))
  str(tf.idf)
  tf.idf@x[1:4]
  tf.idf@p[1:5]
  tf.idf@x <- tf.idf@x * rep(idf.wei, diff(tf.idf@p))
  
  marker.list <- lapply(marker.list, function(m) {
    intersect(m, rownames(obj))}) 
  marker.list <- marker.list[lapply(marker.list, length)>0]
  # Check whether the markers in the seurat object:
  #missed.markers <- setdiff(unlist(marker.list), rownames(obj))
  #print(c("The following markers are missing from the dataset", missed.markers, 
  #  "please make sure each celltype has at least one marker in the dataset!"))

  tf.idf <- as.matrix(tf.idf)
  ti <- lapply(marker.list, function(markers){
    rowSums(tf.idf[, markers, drop=FALSE])
  }) %>% as.data.frame()

  annotation <- c() 
  for (i in 1:nrow(ti)) {
    ind <- which.max(ti[i, ])
    annotation <- c(annotation, colnames(ti)[ind])
  }
  obj$annotation_cell <- setNames(annotation, colnames(obj))

  # Check the correspondence between the annotation and the Louvain clusters
  tab <- plyr::count(obj@meta.data, vars=c(classification, "annotation_cell"))
  tab[, "classification"] <- tab[, classification]
  tab %>% group_by(classification) %>% top_n(n=1, wt="freq")
  top1 <- tab %>% group_by(classification) %>% top_n(n=1, wt=freq)
  top3 <- tab %>% group_by(classification) %>% top_n(n=3, wt=freq)
  # For tie situation
  clr <- unique(top1[,1][[1]])
  first <- lapply(clr, function(i) {
        print(top1[top1[,1][[1]] %in% i, ][1,])
  })
  top1 <- do.call(rbind, first)

  obj@meta.data[, "annotation_cluster"] <- obj@meta.data[, classification]
  new.cluster.ids <- top1$annotation_cell ####
  obj <- changeLabel(obj, "annotation_cluster", new.cluster.ids)

  return(list(obj=obj, top3ann=top3))
}



#' Make a stacked VlnPlot for that celltypes acquired from the annotation
#' @param obj a seurat object
#' @param marker.list a list of markers for different cell types
#' @param classification classification to be annotated
#' @param plot.margin  plot margin
#' @export stackedVlnPlot
stackedVlnPlot <- function(obj, marker.list, plot.margin = unit(c(0.9, 1, 1, 2), "cm")) {
  Idents(obj) <- obj$annotation_cluster
  markers <- marker.list[levels(obj$annotation_cluster)] %>% unlist
  markers <- intersect(markers, rownames(obj))
  p <- VlnPlot(obj, features=markers, stack = T, sort = F, flip = TRUE) + NoLegend() +
       theme(legend.position = "null", 
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = rel(0)), 
          plot.margin = unit(c(0.9, 1, 1, 2), "cm"))
  return(p)
}


