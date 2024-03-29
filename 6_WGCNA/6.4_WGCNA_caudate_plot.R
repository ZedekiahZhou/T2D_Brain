library(WGCNA)
library(extrafont)
loadfonts()
rm(list=ls())
setwd("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/")
options(stringsAsFactors = FALSE)
tmp <- load("/data/MyProgram/Final_diabrain/7.interdata/6_WGCNA/Brain_Caudate_covar_WGCNAres_v1.0.RData")
bw_label2color <- read.table("/data/MyProgram/Final_diabrain/7.interdata/6_WGCNA/Brain_Caudate_bw_label2color.txt", header = T)
id2symbol <- read.table("/data/MyProgram/Final_diabrain/1.clean/id2symbol.tab", header = T)
mytissue <- "Caudate"
cutoff <- 0.05

ttestpvalue <- function(x, bi_var) {
    tmp <- cbind(vx = x, vbi = bi_var)
    t.test(vx ~ vbi, data = tmp)$p.value
}

#----------------------------------------------------------------------------------------------------------------
# I. module and diabetes
#----------------------------------------------------------------------------------------------------------------
MEs <- moduleEigengenes(datExpr, bwnetColors)$eigengenes
MEs <- MEs[, colnames(MEs)!="MEgrey"]
MEs <- orderMEs(MEs)
bwmoduleTraitCor <- cor(MEs, datTraits, use="p")
bwmoduleTraitPvalue <- corPvalueStudent(bwmoduleTraitCor, nSamples)

# nominal_p <- sapply(MEs, ttestpvalue, datTraits$DIABETES)
bwmoduleTraitPvalue[, "DIABETES"] <- sapply(MEs, ttestpvalue, datTraits$DIABETES)
diab_padj <- p.adjust(bwmoduleTraitPvalue[, "DIABETES"], method = "BH")



textMatrix <- paste(signif(bwmoduleTraitCor, 2), "\n(", signif(bwmoduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(bwmoduleTraitCor)
#sizeGrWindow(7, 6)
labeledHeatmap(Matrix = bwmoduleTraitCor, xLabels = colnames(bwmoduleTraitCor), 
               yLabels = rownames(bwmoduleTraitCor), ySymbols = rownames(bwmoduleTraitCor),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = F, 
               cex.text = 0.5, zlim = c(-1,1))
pdf(paste0("./Module_Trait_", mytissue, ".pdf"), family = "Arial", width = 3, height = 6)
par(mar = c(2, 7, 2, 1.5), ps = 10)
labeledHeatmap(Matrix = bwmoduleTraitCor[, "DIABETES", drop = FALSE], xLabels = "T2D", xLabelsAngle = 0, xLabelsAdj = 0.5, 
               yLabels = rownames(bwmoduleTraitCor), ySymbols = rownames(bwmoduleTraitCor), 
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix[, 6], setStdMargins = F, 
               cex.text = 1.0, cex.lab = 1.0, zlim = c(-1,1))
dev.off()
# embed_fonts(paste0("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/Module_Trait_", mytissue, ".pdf"))
cor_mods <- sub("^ME", "", rownames(bwmoduleTraitPvalue[bwmoduleTraitPvalue[, "DIABETES"] < cutoff, ]))

# generate a summary table
all(colnames(bwmoduleTraitCor) == colnames(bwmoduleTraitPvalue))
diab_module_cor <- data.frame(row.names(bwmoduleTraitCor), bwmoduleTraitCor[, "DIABETES"], bwmoduleTraitPvalue[, "DIABETES"], diab_padj)
colnames(diab_module_cor) <- c("Module", "Correlation", "Nominal P-value", "Adjust P-value")
module_number <- as.data.frame(table(bwnetColors))
colnames(module_number) <- c("Module", "Size")
module_number$Module <- paste0("ME", module_number$Module)
diab_module_cor <- merge(diab_module_cor, module_number)

write.table(diab_module_cor, paste0("./Module_Trait_", mytissue, ".txt"), quote = F, sep = "\t", row.names = F)

# cor_mods <- cor_mods[cor_mods != "grey"]

#-------------------------------------------------------------------------------------------------------------
# II. Gene Significance and Moudule Membership
#-------------------------------------------------------------------------------------------------------------
diab <- as.data.frame(datTraits$DIABETES)
names(diab) <- "DIABETES"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
# MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
# names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, diab, use = "p"))
# GSPvalue <- as.data.frame(sapply(datExpr, FUN = ttestpvalue, diab[, 1]))
names(geneTraitSignificance) <- paste("GS.", names(diab), sep = "")
# names(GSPvalue) <- paste("p.GS.", names(diab), sep = "")

pdf("./GSvsMM_all_mod.pdf", family = "Arial", width = 9, height = 9)
for (module in modNames) {
    column <- match(module, modNames)
    moduleGenes <- bwnetColors==module
    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for T2D",
                       main = paste("Gene significance vs. module membership\n"),
                       abline = TRUE, abline.color = module, 
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/GSvsMM_all_mod.pdf")

# cor_mods <- names(diab_padj)[diab_padj < 0.05]
# cor_mods <- gsub("ME", "", cor_mods)
pdf("./GSvsMM_interest_mod.pdf", family = "Arial", width = 12, height = 8)
# split the canvas and prepare a out margin for title
layout(matrix(c(1:6), 2, 3, byrow = T))
par(mar = c(5, 5, 4, 2), oma = c(0,1,3,0), ps = 16)
for (module in cor_mods) {
    column <- match(module, modNames)
    moduleGenes <- bwnetColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("MM in", module, "module"),
                       ylab = "GS for T2D",
                       #main = paste("Module membership vs. gene significance\n"),
                       abline = TRUE, abline.color = module, 
                       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = module)
}
mtext(expression(bold("Gene significance vs. module membership")), side = 3, line = 0, outer = T, cex = 1)
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/GSvsMM_interest_mod.pdf")

# generate a summary table
mm <- cbind(Gene = rownames(geneModuleMembership), geneModuleMembership)
gs <- cbind(Gene = rownames(geneTraitSignificance), geneTraitSignificance)
gm <- data.frame(Gene = colnames(datExpr), Module = bwnetColors)
gm_and_mm <- merge(gm, gs, by = "Gene" )
gm_and_mm <- merge(gm_and_mm, mm, by = "Gene")
gm_and_mm <- gm_and_mm[gm_and_mm$Module != "grey", ]
gm_and_mm <- gm_and_mm[order(gm_and_mm$Module), ]
gm_and_mm$Gene <- id2symbol$Description[match(gm_and_mm$Gene, id2symbol$Name)]
write.table(gm_and_mm, "./GS_and_MM.txt", sep = "\t", quote = F, row.names = F)

#-------------------------------------------------------------------------------------------------------------
# III. plot GS across module
#-------------------------------------------------------------------------------------------------------------
pdf("./GS_across_module.pdf", family = "Arial", width = 6, height = 6)
par(mar = c(2, 3, 2, 2), ps = 10)
plotModuleSignificance(geneTraitSignificance[[1]], bwnetColors, xaxt="n")
# labs <- paste(names(table(mtcars$cyl)), "cylinders")
# labs <- names(table(bwnetColors))
# text(cex=1, x=c(1:15)-.25, y=-1.25, labs, xpd=TRUE, srt=45)
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/GS_across_module.pdf")
