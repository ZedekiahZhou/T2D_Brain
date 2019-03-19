
# ====================================================================================================
# Author:   Aaron
# Function: WGCNA for caudate
#           Several parameters of blockwiseModules can be changed to get better performance
#           1. high minKMEtostay will make module smaller and module genes with higher KME in each module, lowKMEtostay may make modult split
#           2. unsign network will make smaller modules, genes in singed network corralated with status in the same direction
#           3. bicor is more robust than pearson, but bicor may not suit, choose pearson
# Version:  1.0
# Date:     Oct 15, 2018
# ====================================================================================================

rm(list=ls())
library(DESeq2)
library(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = F)

# Use r-log transformed data from DESeq2, covariates such as SEX, AGE, BMI, RIN were also regress out by a linear regression
tmethod <- "covar_prepare"    
tpath <- paste0("/data/MyProgram/Final_diabrain/7.scripts/6_WGCNA/", tmethod)
id2symbol <- read.table("/data/MyProgram/Final_diabrain/1.clean/id2symbol.tab", stringsAsFactors = F, header = T)
setwd(tpath)

#-------------------------------------------------------------------------------------------------
# I. Check data
#-------------------------------------------------------------------------------------------------
# specify which tissue to research
mytissue <- "Brain_Caudate"
fname <- dir("./")
tname <- sub("\\.RData", "", fname)

#for (mytissue in tname) {
    load(paste0(mytissue, ".RData"))
    datTraits <- mycol[, c(-1, -2, -9)]
    datTraits$SEX <- as.numeric(as.character(datTraits$SEX))
    datTraits$RACE <- as.numeric((as.character(datTraits$RACE)))
    datTraits$DIABETES <- as.numeric(as.character(datTraits$DIABETES))
    rownames(datTraits) <- make.names(mycol$SAMPID)
    datExpr = as.data.frame(t(dds_residual))
    all(rownames(datExpr)==rownames(datTraits))
    
    nGenes <- dim(datExpr)[2]
    nSamples <- dim(datExpr)[1]
    
    #-------------------------------------------------------------------------------------------------
    # II. Pick Soft Threshold
    #-------------------------------------------------------------------------------------------------
    # Choose power
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    # Call the network topology analysis function
    system.time(sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize = 30000, RsquaredCut = 0.90))
    # Plot the results:
    #sizeGrWindow(9, 5)
    pdf(file = paste0("../../../4.plots/WGCNA_Caudate/softpower_",mytissue,".pdf"), width = 9, height = 5);
    par(mfrow = c(1,2));
    cex1 = 0.9;
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity "))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
    
    #-------------------------------------------------------------------------------------------------
    # III. Net Construction
    #-------------------------------------------------------------------------------------------------
    collectGarbage();
    spower <- sft$powerEstimate
    if (is.na(spower)) {
        spower <- sft$fitIndices[sft$fitIndices[,2]==max(sft$fitIndices[,2]), 1]
    }
    
    aconnect <- softConnectivity(datExpr, type = "signed", power = spower, blockSize = 30000)
    # One-step
    bwnet <- blockwiseModules(datExpr, maxBlockSize = 30000,
                              power = spower, TOMType = "unsigned",
                              minModuleSize = 40, reassignThreshold = 0,
                              mergeCutHeight = 0.15, numericLabels = TRUE,
                              corType = "pearson", networkType = "signed", 
                              minKMEtoStay = 0.7, 
                              saveTOMs = FALSE, verbose = 3)
    table(bwnet$colors)
   
    bwnetLabels <- bwnet$colors
    bwnetColors <- labels2colors(bwnetLabels)
    bw_label2color <- data.frame(Labels = bwnetLabels,  Colors = bwnetColors)
    bw_label2color <- unique(bw_label2color)
    write.table(bw_label2color, file = paste0("../", mytissue, "_bw_label2color.txt"), quote = F, sep = "\t", row.names = F)
    
    pdf(file = paste0("../../../4.plots/WGCNA_Caudate/DendroColor_",mytissue,".pdf"), width = 12, height = 8);
    plotDendroAndColors(bwnet$dendrograms[[1]], bwnetColors,
                        sub("Brain_", "", mytissue), dendroLabels = FALSE, hang = 0.03,
                        main = paste0("Gene dendrogram and module colors of ", sub("Brain_", "", mytissue)),
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
    mod2gene <- data.frame(ModID = bwnetColors, geneID = colnames(datExpr), stringsAsFactors = F)
    mod2gene <- mod2gene[order(mod2gene$ModID), ]
    mod2gene <- data.frame(Mod = paste0("ME", mod2gene$ModID), mod2gene, stringsAsFactors = F)
    mod2gene <- mod2gene[, -2]
    
    tmp <- dir("../../../2.results/Cov_optimal/result05/res/", full.names = T)
    DEgene <- read.table(tmp[grep(mytissue, tmp)], header = T, stringsAsFactors = F)
    togene <- data.frame(Mod = 'TotalGene', geneID = DEgene$Gene)
    upgene <- data.frame(Mod = "UpGene", geneID = DEgene$Gene[DEgene$log2FoldChange > 0])
    downgene <- data.frame(Mod = "DownGene", geneID = DEgene$Gene[DEgene$log2FoldChange < 0])
    mod2gene <- rbind(mod2gene, togene, upgene, downgene)
    
    
    mod2gene$geneID <- id2symbol$Description[match(mod2gene$geneID, id2symbol$Name)]
    
    save(mod2gene, bwnet, datExpr, datTraits, bwnetLabels, bwnetColors, nGenes, nSamples, bw_label2color, file = paste0("../", mytissue, "_covar_WGCNAres_v1.0.RData"))
#}
