
# CellAnnotatoRSeurat
**Annotate cells/cluster based on a marker list**

# Quick guide

### 1. Load the package:
**Bash**

    cd ~/bin/
    git clone https://github.com/Maj18/CellAnnotatoRSeurat.git

**R**

    devtools::load_all("~/bin/CellAnnotatoRSeurat")
    #devtools::install_github("Maj18/CellAnnotatoRSeurat")
    library(Matrix)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    
### 2. Import data

    DIR <- "~/bin/CellAnnotatoRSeurat/example/"
    INDIR <- paste0(DIR, "/data/")
    obj <- readRDS(paste0(INDIR, "gw20_integrated.RDS"))
    
### 3. Define a marker list

    marker.list <- list(NE=c("VIM", "HMGA2", "ARHGAP28", "ALCAM"), NP=c("NES", "ASCL1", "NBPF10", "NBPF15", "PLAGL1", "SFTA3"),
                    Neuron=c("STMN2", "SYT1", "SLC32A1", "GAD2", "HDC", "SLC17A6"), OPC=c("PDGFRA", "GPR75-ASB3"),
                    OL=c("OLIG2", "GPR75-ASB3", "MBP", "MOG"), Astrocyte=c("AQP4", "AGT", "FAM107A"), 
                    Ependy=c("CCDC153", "CCDC74A", "CCDC74B", "LRTOMT"), Microglia=c("AIF1", "CX3CR1","ITGAM", "PTPRC"), 
                    Endoth=c("CLDN5", "FLT1", "SLC38A5"), Mural=c("NDUFA4L2", "PDGFRB"), VLMC=c("PTGDS", "COL1A1"))

### 4. Annotation by cell and by cluster, based on provided markers above

classification="integrated_snn_res.0.1"
rslt <- annOnMarker(obj=obj, marker.list, classification=classification)
obj <- rslt$obj

### 5. Check results

**Plot annotation by cell, by cluster, as well as by the original classification side by side**

    cowplot::plot_grid(
        DimPlot(obj, group.by="annotation_cell", label=T, label.box=T, repel=T) + NoLegend(),
        DimPlot(obj, group.by=classification, label=T, label.box=T, repel=T) + NoLegend(),
        DimPlot(obj, group.by="annotation_cluster", label=T, label.box=T, repel=T) + NoLegend(),
        nrow=2
        )
        
**To check the quality of the annotation by cluster, we can make a stacked VlnPlot** 

    stackedVlnPlot(obj, marker.list)
    
**To check other alternative top annotation for each cluster**

    rslt$top3ann

