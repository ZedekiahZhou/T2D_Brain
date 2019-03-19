# ====================================================================================================
# Author:   Aaron
# Function: randomly choose subset of T2D samples with matched control samples, use DESeq2 to identify
#           diferentially expressed genes. 
# Version:  1.0
# Date:     Jun 3, 2019
# ====================================================================================================

#Note: The ComplexHeatmap package uses grid graphics, and it is often the case that you need to (print) your (Heatmap) output for it to show up in the pdf.
rm(list=ls())
library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
fpath <- "/data/MyProgram/Final_diabrain/7.scripts/6_WGCNA/covar_prepare/"
setwd(fpath)

#set the initial value
#dis <- "euclidean"

cutoff <- 0.05
tmp <- load("Brain_Caudate.RData")
tname <- "Caudate"

# tmp <- load("Brain_Hippocampus.RData")
# tname <- "Hippocampus"

    print(paste(tname, "is begin at", Sys.time()))
    saminfo_file <- dir("/data/MyProgram/Final_diabrain/1.clean/matched_sam_optimal/", full.names = T)
    tsaminfo <- read.table(saminfo_file[grep(tname, saminfo_file)], 
                           header = T, sep="\t", quote="", stringsAsFactors = F) 
    colnames(tsaminfo) <- c("SUBJID", "SAMPID", "SEX", "AGE", "RACE", "BMI", "RIN", "DIABETES")
    anno_info <- tsaminfo
    anno_info$SEX <- as.factor(ifelse(tsaminfo$SEX==1, "Male", "Female"))
    anno_info$RACE <- as.factor(ifelse(tsaminfo$RACE==2, "Black", "White"))
    anno_info$DIABETES <- as.factor(ifelse(tsaminfo$DIABETES==0, "Control", "T2D"))
    
    #Heatmap--------------------------------------------------------------------
    
    #Check if the rownames of tsaminfo and colnames of tdata are the same
    tres <- results(tdds, alpha=cutoff)
    DEix <- (!is.na(tres$padj))&(tres$padj<cutoff)
    tdata <- t(scale(t(dds_residual[DEix, ])))
    #tdata <- tdata[, T2D]
    rownames(tsaminfo) <- make.names(tsaminfo$SAMPID)
    ix <- all(rownames(tsaminfo) == colnames(tdata))
    ix
    if (!ix) {tdata <- tdata[, rownames(tsaminfo)]}
    
    
    
    #No main tile, so use the column title to define the main tile
    #clus <- "complete"
    myCols <- brewer.pal(n=3, name="RdBu")
    annoCol <- brewer.pal(n= 12, name = "Paired")
    annot_data <- anno_info[, c("DIABETES", "AGE", "SEX", "BMI", "RACE", "RIN")]
    ha <- HeatmapAnnotation(annot_data, col = list(DIABETES = c("Control" = "white", "T2D"=annoCol[1]), 
                                                   AGE = colorRamp2(c(min(annot_data$AGE), max(annot_data$AGE)), c("white", annoCol[2])),
                                                   SEX = c("Male" = annoCol[3], "Female" = annoCol[4]),
                                                   BMI = colorRamp2(c(min(annot_data$BMI), max(annot_data$BMI)), c("white", annoCol[5])), 
                                                   RACE = c("White" = annoCol[7], "Black" = annoCol[8]),
                                                   RIN = colorRamp2(c(min(annot_data$RIN), max(annot_data$RIN)), c("white", annoCol[9]))
                                                  ))
    cskip <- 0.01
    cmin <- quantile(as.vector(tdata), probs=cskip)
    cmax <- quantile(as.vector(tdata), probs=1-cskip)
    pdf(paste0("../../../4.plots/heatmap/", tname, ".pdf"), width=9, height = 9)
    #dev.new()
    a <- Heatmap(tdata, col = colorRamp2(c(cmin, 0, cmax), myCols), 
            column_title =paste("Type 2 Diabetes-associated gene expression in", tname),
            column_title_side = "bottom", row_title = "DE Genes", 
            show_column_names = FALSE, 
            show_row_names = FALSE, 
            clustering_distance_columns = "pearson", 
            clustering_distance_rows = "euclidean", 
            clustering_method_rows = "ward.D2", 
            clustering_method_columns = "ward.D",
            column_dend_height = unit(2, "cm"),
            row_dend_width = unit(2, "cm"),
            heatmap_legend_param = list(title="expression"), 
            top_annotation = ha)
    print(a)
    
    # Use fisher test to test if there is difference in the diabetes status of the two branches of the Heatmap sample cluster result
    tcolden <- column_dend(a)
    ## get the sample id in the order of dendrogram
    n <- attributes(tcolden)$members
    ## the cut point 
    cutpoint <- attributes(tcolden[[1]])$members
    ## sample id of left branches 
    tT2D_id <- order.dendrogram(tcolden)[1:cutpoint]
    ## sample id of right branches 
    tCtrl_id <- order.dendrogram(tcolden)[cutpoint+1:n]
    ## get the diabetes status and perform fisher test
    tT2D_diab <- table(anno_info$DIABETES[tT2D_id])
    tCtrl_diab <- table(anno_info$DIABETES[tCtrl_id])
    p_value <- fisher.test(cbind(tT2D_diab, tCtrl_diab))$p.value
    plot.new()
    text(0.5, 0, paste0("Fisher test P-value: ", format(p_value, digits = 3, scientific = TRUE)))
    
    dev.off()


