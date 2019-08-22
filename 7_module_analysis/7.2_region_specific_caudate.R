
# ====================================================================================================
# Author:   Aaron
# Function: region specific test for caudate
# Version:  1.0
# Date:     Mar 17, 2019
# ====================================================================================================
rm(list=ls())
library(pSI)
library(pSI.data)
library(WGCNA)
library(extrafont)
loadfonts()
options(stringsAsFactors = FALSE)
setwd("/data/MyProgram/Final_diabrain/7.interdata/7_GWAS_REGION/script/")
tmp <- load("/data/MyProgram/Final_diabrain/7.interdata/6_WGCNA/Brain_Caudate_covar_WGCNAres_v1.0.RData")
id2symbol <- read.table("/data/MyProgram/Final_diabrain/1.clean/id2symbol.tab", stringsAsFactors = F, header = T)
celltype_name_transform <- read.table("./celltype_name_transform.tab", header = T, sep = "\t")


#----------------------------------------------------------------------------------------------------------------
# I. region enrichment
#----------------------------------------------------------------------------------------------------------------
modgene_symbol <- split(mod2gene$geneID, mod2gene$Mod )

data(human);data(mouse)

region_pSI <- list()
celltype_pSI <- list()

for (i in 1:length(modgene_symbol)) {
    region_pSI[[i]] <- fisher.iteration(human$young.adulthood$psi.out, modgene_symbol[[i]], p.adjust = F)
    celltype_pSI[[i]] <- fisher.iteration(mouse$psi.out, modgene_symbol[[i]], background = "mouse.human", p.adjust = F)
}
names(region_pSI) <- names(modgene_symbol)
names(celltype_pSI) <- names(modgene_symbol)
# filter out regions in mouse cell type data
celltype_pSI <- lapply(celltype_pSI, function(x) {x[match(celltype_name_transform$ROW_NAME, rownames(x)), ]})

# convert list result from fisher.interation to data.frame
list2df <- function(x, idx) {
    tmp <- x[, idx]
    names(tmp) <- rownames(x)
    tmp
}

# tidy up for out put
nominal_plist <- list(r0001 = sapply(region_pSI, list2df, 4),
                      r001 = sapply(region_pSI, list2df, 3),
                      r01 = sapply(region_pSI, list2df, 2),
                      c0001 = sapply(celltype_pSI, list2df, 4),
                      c001 = sapply(celltype_pSI, list2df, 3),
                      c01 = sapply(celltype_pSI, list2df, 2))

# bonferroni adjust of p value
my_bonf_adj <- function(x) {
    n <- dim(x)[1]*dim(x)[2]
    res <- x*n
    res[res>1] <- 1
    rownames(res) <- paste0(rownames(x), "_adj")
    res <- rbind(x, res)
    res <- res[order(rownames(res)), ]
    res <- as.data.frame(t(res))
    res
}

bonf_plist <- list(bonf_r0001 = my_bonf_adj(nominal_plist$r0001),
                   bonf_r001 = my_bonf_adj(nominal_plist$r001),
                   bonf_r01 = my_bonf_adj(nominal_plist$r01),
                   bonf_c0001 = my_bonf_adj(nominal_plist$c0001),
                   bonf_c001 = my_bonf_adj(nominal_plist$c001),
                   bonf_c01 = my_bonf_adj(nominal_plist$c01))
write.table(rbind(bonf_plist$bonf_r0001, bonf_plist$bonf_r001, bonf_plist$bonf_r01), 
            file = "/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/region_specificity.tab", quote = F, sep = "\t")
write.table(rbind(bonf_plist$bonf_c0001, bonf_plist$bonf_c001, bonf_plist$bonf_c01), 
            file = "/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/celltype_specificity.tab", quote = F, sep = "\t")


#----------------------------------------------------------------------------------------------------------------
# I. plot the enrichment heatmap
#----------------------------------------------------------------------------------------------------------------
# idx <- c(colnames(r0001)[grep("ME", colnames(r0001))], "UpGene", "DownGene")
# region_plot <- -log10(r0001)[, c(2:14, 16, 1)]
idx <- c(rownames(bonf_plist$bonf_r0001)[grep("ME", rownames(bonf_plist$bonf_r0001))])
idx <- idx[idx != "MEgrey"]
idy <- c(colnames(bonf_plist$bonf_r0001)[grep("_adj", colnames(bonf_plist$bonf_r0001))])
region_plot <- -log10(bonf_plist$bonf_r0001)[idx, idy]
rownames(region_plot) <- sub("ME", "", rownames(region_plot))
text_region <- matrix("", nrow = dim(region_plot)[1], ncol = dim(region_plot)[2])
for (i in 1:dim(region_plot)[1]) {
    for (j in 1:dim(region_plot)[2]) {
        if (region_plot[i, j] > 4) {
            text_region[i, j] <- "****"
        } else if (region_plot[i, j] > 3) {
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
idy <- c(colnames(bonf_plist$bonf_c0001)[grep("_adj", colnames(bonf_plist$bonf_c0001))])
celltype_plot <- -log10(bonf_plist$bonf_c0001)[idx, idy]
text_celltype <- matrix("", nrow = dim(celltype_plot)[1], ncol = dim(celltype_plot)[2])
for (i in 1:dim(celltype_plot)[1]) {
    for (j in 1:dim(celltype_plot)[2]) {
        if (celltype_plot[i, j] > 4) {
            text_celltype[i, j] <- "****"
        } else if (celltype_plot[i, j] > 3) {
            text_celltype[i, j] <- "***"
        } else if (celltype_plot[i, j] > 2) {
            text_celltype[i, j] <- "**"
        } else if (celltype_plot[i, j] > -log10(0.05)) {
            text_celltype[i, j] <- "*"
        }
    }
}
max_celltype <- max(celltype_plot)

# pdf("../../../4.plots/WGCNA_Caudate/region_specificity.pdf", family ="Arial", height = 12, width = 12)
# # vertical 
# layout(matrix(c(rep(1,3), 2), ncol = 1))
# #layout.show(2)
# cell_ylab <- celltype_name_transform$Description[match(sub("_adj", "", colnames(celltype_plot)), celltype_name_transform$ROW_NAME)]
# par(mar = c(0, 60, 3, 3), ps=7)
# labeledHeatmap(Matrix = t(celltype_plot), xLabels = rep("", dim(celltype_plot)[1]), #xLabelsPosition = "top", #xSymbols = colnames(celltype_plot), 
#                yLabels = cell_ylab,  textMatrix = t(text_celltype),
#                colorLabels = FALSE, colors = blueWhiteRed(50)[26:50],  setStdMargins = F, 
#                cex.text = 1.2, cex.lab = 1.5, zlim = c(0, max_celltype))
# par(mar = c(9, 60, 0, 3))
# labeledHeatmap(Matrix = t(region_plot), xLabels = rownames(region_plot), xSymbols = rownames(region_plot),
#                yLabels = sub("\\.[[:print:]]+", "", colnames(region_plot)), textMatrix = t(text_region),
#                colorLabels = FALSE, colors = blueWhiteRed(50)[26:50],  setStdMargins = F, 
#                cex.text = 1.2, cex.lab = 1.5, zlim = c(0, max_region))
# 
# # horizonical 
# layout(matrix(1, ncol = 1))
# par(mar = c(40, 10, 3, 3))
# merge_plot <- cbind(region_plot, celltype_plot)
# ylabel <- rownames(region_plot)
# xlabel <- c(sub("\\.[[:print:]]+", "", colnames(region_plot)), cell_ylab)
# text_merge <- cbind(text_region, text_celltype)
# labeledHeatmap(Matrix = merge_plot, xLabels = xlabel, yLabels = ylabel, ySymbols = ylabel,textMatrix = text_merge,
#                colorLabels = FALSE, colors = blueWhiteRed(50)[26:50],  setStdMargins = F, 
#                cex.text = 1.1, cex.lab = 1.3, zlim = c(0, max_celltype), verticalSeparator.x = 6, verticalSeparator.ext = 0)
# dev.off()

pdf("../../../4.plots/WGCNA_Caudate/region_specificity_new.pdf", family ="Arial", height = 12, width = 8)
# vertical
layout(matrix(c(rep(1,3), 2), ncol = 1))
#layout.show(2)
cell_ylab <- celltype_name_transform$Description[match(sub("_adj", "", colnames(celltype_plot)), celltype_name_transform$ROW_NAME)]
par(mar = c(0, 30, 3, 3), ps=15)
labeledHeatmap(Matrix = t(celltype_plot), xLabels = rep("", dim(celltype_plot)[1]), #xLabelsPosition = "top", #xSymbols = colnames(celltype_plot),
               yLabels = cell_ylab,  textMatrix = t(text_celltype),
               colorLabels = FALSE, colors = blueWhiteRed(50)[26:50],  setStdMargins = F,
               cex.text = 1.2, cex.lab = 1.5, zlim = c(0, max_celltype))
par(mar = c(9, 30, 0, 3))
labeledHeatmap(Matrix = t(region_plot), xLabels = rownames(region_plot), xSymbols = rownames(region_plot),
               yLabels = sub("\\.[[:print:]]+", "", colnames(region_plot)), textMatrix = t(text_region),
               colorLabels = FALSE, colors = blueWhiteRed(50)[26:50],  setStdMargins = F,
               cex.text = 1.2, cex.lab = 1.5, zlim = c(0, max_region))
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/region_specificity.pdf")
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
save(modgene_symbol, region_pSI, celltype_pSI, bonf_plist, nominal_plist, valuable_subMod001, valuable_subMod01, file = "./caudate/region_v1.0.RData")
