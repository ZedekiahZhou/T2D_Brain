
# ====================================================================================================
# Author:   Aaron
# Function: region specific test for Hippocampus
# Version:  1.0
# Date:     Mar 17, 2019
# ====================================================================================================
rm(list=ls())
library(pSI)
library(pSI.data)
library(WGCNA)
setwd("/data/MyProgram/Final_diabrain/7.scripts/7_GWAS_REGION/script/")
tmp <- load("/data/MyProgram/Final_diabrain/7.scripts/6_WGCNA/Brain_Hippocampus_covar_WGCNAres_v1.0.RData")
id2symbol <- read.table("/data/MyProgram/Final_diabrain/1.clean/id2symbol.tab", stringsAsFactors = F, header = T)
options(stringsAsFactors = FALSE)

#----------------------------------------------------------------------------------------------------------------
# I. region enrichment
#----------------------------------------------------------------------------------------------------------------
modgene_symbol <- split(mod2gene$geneID, mod2gene$Mod )

data(human);data(mouse)

region_pSI <- list()
celltype_pSI <- list()

for (i in 1:length(modgene_symbol)) {
    region_pSI[[i]] <- fisher.iteration(human$young.adulthood$psi.out, modgene_symbol[[i]])
    celltype_pSI[[i]] <- fisher.iteration(mouse$psi.out, modgene_symbol[[i]], background = "mouse.human")
}
names(region_pSI) <- names(modgene_symbol)
names(celltype_pSI) <- names(modgene_symbol)

list2df <- function(x, idx) {
    tmp <- x[, idx]
    names(tmp) <- rownames(x)
    tmp
}

r0001 <- sapply(region_pSI, list2df, 4)
r001 <- sapply(region_pSI, list2df, 3)
r01 <- sapply(region_pSI, list2df, 2)
r05 <- sapply(region_pSI, list2df, 1)
c0001 <- sapply(celltype_pSI, list2df, 4)
c001 <- sapply(celltype_pSI, list2df, 3)
c01 <- sapply(celltype_pSI, list2df, 2)
c05 <- sapply(celltype_pSI, list2df, 1)


#----------------------------------------------------------------------------------------------------------------
# I. plot the enrichment heatmap
#----------------------------------------------------------------------------------------------------------------
idx <- c(colnames(r0001)[grep("ME", colnames(r0001))])
idx <- idx[idx != "MEgrey"]
# region_plot <- -log10(r0001)[, c(2:14, 16, 1)]
region_plot <- -log10(r0001)[, idx]
text_region <- matrix("", nrow = dim(region_plot)[1], ncol = dim(region_plot)[2])
for (i in 1:dim(region_plot)[1]) {
    for (j in 1:dim(region_plot)[2]) {
        if (region_plot[i,j] > 3) {
            text_region[i, j] <- "***"
        } else if (region_plot[i,j] > 2) {
            text_region[i, j] <- "**"
        } else if (region_plot[i,j] > -log10(0.05)) {
            text_region[i,j] <- "*"
        }
        
    }
}
max_region <- max(region_plot)

# celltype_plot <- -log10(c0001)[, c(2:14, 16, 1)]
celltype_plot <- -log10(c0001)[, idx]
text_celltype <- matrix("", nrow = dim(celltype_plot)[1], ncol = dim(celltype_plot)[2])
for (i in 1:dim(celltype_plot)[1]) {
    for (j in 1:dim(celltype_plot)[2]) {
        if (celltype_plot[i,j] > 3) {
            text_celltype[i,j] <- "***"
        } else if (celltype_plot[i,j] > 2) {
            text_celltype[i,j] <- "**"
        } else if (celltype_plot[i,j] > -log10(0.05)) {
            text_celltype[i,j] <- "*"
        }
    }
}
max_celltype <- max(celltype_plot)
celltype_name_transform <- read.table("./celltype_name_transform.tab", header = T, sep = "\t")

pdf("../../../4.plots/WGCNA_Hippocampus/region_specificity.pdf", height = 12, width = 12)
layout(matrix(c(rep(1,4), 2), ncol = 1))
#layout.show(2)
cell_ylab <- celltype_name_transform$Description[match(rownames(celltype_plot), celltype_name_transform$ROW_NAME)]
par(mar = c(0, 35, 3, 3))
labeledHeatmap(Matrix = celltype_plot, xLabels = rep("", dim(celltype_plot)[2]), #xLabelsPosition = "top", #xSymbols = colnames(celltype_plot), 
               yLabels = cell_ylab,  textMatrix = text_celltype,
               colorLabels = FALSE, colors = blueWhiteRed(50)[26:50],  setStdMargins = F, 
               cex.text = 1.2, cex.lab = 1.2, zlim = c(0, max_celltype))
par(mar = c(9, 35, 0, 3))
labeledHeatmap(Matrix = region_plot, xLabels = colnames(region_plot), xSymbols = colnames(region_plot),
               yLabels = sub("\\.[[:print:]]+", "", rownames(region_plot)), textMatrix = text_region,
               colorLabels = FALSE, colors = blueWhiteRed(50)[26:50],  setStdMargins = F, 
               cex.text = 1.2, cex.lab = 1.2, zlim = c(0, max_region))
dev.off()

#----------------------------------------------------------------------------------------------------------------
# III. extract all submod genes 
#----------------------------------------------------------------------------------------------------------------
# pSI = 0.001
valuable_subMod001 <- data.frame()
for (i in 1:length(modgene_symbol)) {
    tmpModName <- names(modgene_symbol)[i]
    MEhuman <- candidate.overlap(human$young.adulthood$psi.out, modgene_symbol[[i]])$pSi_0.001
    MEmouse <- candidate.overlap(mouse$psi.out, modgene_symbol[[i]])$pSi_0.001
    for (j in 1:ncol(MEhuman)) {
        if (sum(!is.na(MEhuman[, j]))>10) {
            submodName <- paste(tmpModName, colnames(MEhuman)[j], sep = "_")
            submodgene <- MEhuman[, j]
            submodgene <- submodgene[!is.na(submodgene)]
            l <- length(submodgene)
            valuable_subMod001 <- rbind(valuable_subMod001, data.frame(subMod = rep(submodName, l), Gene = submodgene))
        }
    }
    for (j in 1:ncol(MEmouse)) {
        if (sum(!is.na(MEmouse[, j]))>10) {
            submodName <- paste(tmpModName, colnames(MEmouse)[j], sep = "_")
            submodgene <- MEmouse[, j]
            submodgene <- submodgene[!is.na(submodgene)]
            l <- length(submodgene)
            valuable_subMod001 <- rbind(valuable_subMod001, data.frame(subMod = rep(submodName, l), Gene = submodgene))
        }
    }
}

#pSI = 0.01
valuable_subMod01 <- data.frame()
for (i in 1:length(modgene_symbol)) {
    tmpModName <- names(modgene_symbol)[i]
    MEhuman <- candidate.overlap(human$young.adulthood$psi.out, modgene_symbol[[i]])$pSi_0.01
    MEmouse <- candidate.overlap(mouse$psi.out, modgene_symbol[[i]])$pSi_0.01
    for (j in 1:ncol(MEhuman)) {
        if (sum(!is.na(MEhuman[, j]))>10) {
            submodName <- paste(tmpModName, colnames(MEhuman)[j], sep = "_")
            submodgene <- MEhuman[, j]
            submodgene <- submodgene[!is.na(submodgene)]
            l <- length(submodgene)
            valuable_subMod01 <- rbind(valuable_subMod01, data.frame(subMod = rep(submodName, l), Gene = submodgene))
        }
    }
    for (j in 1:ncol(MEmouse)) {
        if (sum(!is.na(MEmouse[, j]))>10) {
            submodName <- paste(tmpModName, colnames(MEmouse)[j], sep = "_")
            submodgene <- MEmouse[, j]
            submodgene <- submodgene[!is.na(submodgene)]
            l <- length(submodgene)
            valuable_subMod01 <- rbind(valuable_subMod01, data.frame(subMod = rep(submodName, l), Gene = submodgene))
        }
    }
}
save(modgene_symbol, region_pSI, celltype_pSI, r0001, c0001, valuable_subMod001, valuable_subMod01, file = "./hippocampus/region_v1.0.RData")
