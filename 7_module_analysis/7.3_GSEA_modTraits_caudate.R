# ====================================================================================================
# Author:   Aaron
# Function: Modules correlation with diabetes
#           Use GSEA to idetify modules enrich in T2D top GWAS gene list.
#           pathway enrichment
# Version:  1.0
# Date:     Jan 10, 2019
# ====================================================================================================
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(WGCNA)
rm(list=ls())
setwd("/data/MyProgram/Final_diabrain/7.scripts/7_GWAS_REGION/script/")
options(stringsAsFactors = FALSE)

# load diabetes and height GWAS gene list
diab_snp_gene <- read.table("2_DIABETES_SNP_GENEs_for_GSEA.tab", stringsAsFactors = F, header = T)
height_snp_gene <- read.table("2_HEIGHT_SNP_GENEs_for_GSEA.tab", stringsAsFactors = F, header = T)
alzheimer_snp_gene <- read.table("2_ALZHEIMER_SNP_GENES_for_GSEA.tab", header = T)
id2symbol <- read.table("/data/MyProgram/Final_diabrain/1.clean/id2symbol.tab", stringsAsFactors = F, header = T)
tmp <- load("/data/MyProgram/Final_diabrain/7.scripts/6_WGCNA/Brain_Caudate_covar_WGCNAres_v1.0.RData")
tmp <- load("/data/MyProgram/Final_diabrain/7.scripts/7_GWAS_REGION/script/caudate/region_v1.0.RData")

#----------------------------------------------------------------------------------------------------------------
# I. module and diabetes
#----------------------------------------------------------------------------------------------------------------
MEs <- moduleEigengenes(datExpr, bwnetColors)$eigengenes
MEs <- orderMEs(MEs)
bwmoduleTraitCor <- cor(MEs, datTraits, use="p")
bwmoduleTraitPvalue <- corPvalueStudent(bwmoduleTraitCor, nSamples)
for (i in 1:ncol(MEs)) {
    tmp <- cbind(me = MEs[, i],  diab = datTraits$DIABETES)
    bwmoduleTraitPvalue[i, "DIABETES"] <- t.test(me ~ diab, data = tmp)$p.value
}

## plot the MEs with variables
textMatrix <- paste(signif(bwmoduleTraitCor, 2), "\n(", signif(bwmoduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(bwmoduleTraitCor)
sizeGrWindow(10, 6)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = bwmoduleTraitCor, xLabels = colnames(bwmoduleTraitCor), 
               yLabels = rownames(bwmoduleTraitCor), ySymbols = rownames(bwmoduleTraitCor),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = F, 
               cex.text = 0.5, zlim = c(-1,1))



#---------------------------------------------------------------------------------------------------------------
# II. Overlap with DEgenes
#---------------------------------------------------------------------------------------------------------------
nMod <- length(table(bwnet$colors))
Modname <- names(table(bwnetColors))
DEgene <- read.table("../../../2.results/Cov_optimal/result05/res/Brain_Caudate.res", header = T, stringsAsFactors = F)
togene <- DEgene$Gene
upgene <- DEgene$Gene[DEgene$log2FoldChange > 0]
downgene <- DEgene$Gene[DEgene$log2FoldChange < 0]

modgene <- split(names(datExpr), f = as.factor(bwnetColors))
updens <- sapply(modgene, FUN = function(x) {length(intersect(upgene, x))/length(x)})
downdens <- sapply(modgene, FUN = function(x) {length(intersect(downgene, x))/length(x)})
upprop <- sapply(modgene, FUN = function(x) {length(intersect(upgene, x))/length(upgene)})
downprop <- sapply(modgene, FUN = function(x) {length(intersect(downgene, x))/length(downgene)})

up_mod <- enricher(id2symbol$Description[match(upgene, id2symbol$Name)], TERM2GENE = mod2gene, maxGSSize = 10000, pvalueCutoff = 1, qvalueCutoff = 1)
down_mod <- enricher(id2symbol$Description[match(downgene, id2symbol$Name)], TERM2GENE = mod2gene, maxGSSize = 10000, pvalueCutoff = 1, qvalueCutoff = 1)
write.table(up_mod@result, "/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/upmod_enrich.tab", sep = "\t", quote = F, row.names = F)
write.table(down_mod@result, "/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/downmod_enrich.tab", sep = "\t", quote = F, row.names = F)
#y <- GSEA(diab_snpgenelist, TERM2GENE = featuretogene,  pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000)

# Summary
all(rownames(bwmoduleTraitCor)==rownames(bwmoduleTraitPvalue))
tmp1 <- data.frame(Name = rownames(bwmoduleTraitCor), diab_cor = bwmoduleTraitCor[, "DIABETES"], 
                  diab_pvalue = bwmoduleTraitPvalue[, "DIABETES"])
all(names(todens)==names(updens)& names(updens)==names(downdens))
tmp2 <- data.frame(Name = paste0("ME", names(todens)), partsofupinmod = upprop, partsofdowninmod = downprop, parts_in_up = updens, part_in_down = downdens)
mod_summary <- merge(tmp1, tmp2, by = "Name")
write.table(mod_summary, "/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/mod_summary.tab", row.names = F, sep = "\t", quote = F)



#--------------------------------------------------------------------------------------------------------------
# III. GSEA for GWAS
#--------------------------------------------------------------------------------------------------------------
## remove duplicate function
rm_dup <- function(snp_gene) {
    tmp <- table(snp_gene$Name)
    dup_gene <- names(tmp)[tmp > 1]
    dup_gene <- data.frame(dup_gene, used = rep("no", length(dup_gene)), stringsAsFactors = F)
    idx <- logical(dim(snp_gene)[1])
    for (i in 1:dim(snp_gene)[1]) {
        j <- match(snp_gene$Name[i], dup_gene$dup_gene)
        if (is.na(j)) {
            idx[i] <- TRUE
        } else {
            if (dup_gene$used[j] == "no") {
                idx[i] <- TRUE
                dup_gene$used[j] <- "yes"
            } else if (dup_gene$used[j] == "yes") {
                idx[i] <- FALSE
            }
        }
    }
    snp_gene <- snp_gene[idx, ]
}

## remove duplicate genes
diab_snp_gene <- rm_dup(diab_snp_gene)
height_snp_gene <- rm_dup(height_snp_gene)
alzheimer_snp_gene <- rm_dup(alzheimer_snp_gene)
## the p.value of top GWAS gene of diabetes is too small for R, so it is calculate as 0 
diab_snp_gene$minus_log10P[1] <- 346.87
alzheimer_snp_gene <- alzheimer_snp_gene[-c(1:7), ]

# prepare genelist for GSEA.P software
write.table(diab_snp_gene[, c("Name", "minus_log10P")], file = "./for_GSEA/diab_genelist.rnk", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(height_snp_gene[, c("Name", "minus_log10P")], file = "./for_GSEA/height_genelist.rnk", sep = "\t", row.names = F, col.names = F, quote = F)

## prepare mod2gene file for GSEA.P software
tmp <- split(mod2gene$geneID, f = mod2gene$Mod)
n <- max(sapply(tmp, FUN = length))
mod_gmx <- matrix(NA, nrow = n+1, ncol = length(tmp))
colnames(mod_gmx) <- names(tmp)
mod_gmx[1, ] <- "module"
for (i in 1:length(tmp)) {
    mod_gmx[2:(length(tmp[[colnames(mod_gmx)[i]]])+1), i] <- tmp[[colnames(mod_gmx)[i]]]
}
mod_gmx <- data.frame(mod_gmx, stringsAsFactors = F)
write.table(mod_gmx, file = "./for_GSEA/module_caudate_v1.0.gmx", quote = F, row.names = F, sep = "\t", na = "")

## prepare submod2gene file for GSEA.P software
valuable_subMod <- valuable_subMod01
tmp <- split(valuable_subMod$Gene, f = valuable_subMod$subMod)
idxtmp <- names(tmp) %in% c("MEturquoise_Cpu.D1_0.01", "MEturquoise_Cpu.D2_0.01", "MEturquoise_Striatum.Young.Adulthood_0.01")
tmp <- tmp[idxtmp]
n <- max(sapply(tmp, FUN = length))
submod_gmx <- matrix(NA, nrow = n+1, ncol = length(tmp))
colnames(submod_gmx) <- names(tmp)
submod_gmx[1, ] <- "sub_module"
for (i in 1:length(tmp)) {
    submod_gmx[2:(length(tmp[[colnames(submod_gmx)[i]]])+1), i] <- tmp[[colnames(submod_gmx)[i]]]
}
submod_gmx <- data.frame(submod_gmx, stringsAsFactors = F)
write.table(submod_gmx, file = "./for_GSEA/sub_module_caudate_v1.0.gmx", quote = F, row.names = F, sep = "\t", na = "")


## Use clusterprofiler to perform GSEA
diab_snpgenelist <- diab_snp_gene$minus_log10P
height_snpgenelist <- height_snp_gene$minus_log10P
alzheimer_snpgenelist <- alzheimer_snp_gene$minus_log10P
names(diab_snpgenelist) <- diab_snp_gene$Name
names(height_snpgenelist) <- height_snp_gene$Name
names(alzheimer_snpgenelist) <- alzheimer_snp_gene$Name

diab_modGSEA <- GSEA(diab_snpgenelist, TERM2GENE = mod2gene, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
diab_modGSEA <- diab_modGSEA@result
height_modGSEA <- GSEA(height_snpgenelist, TERM2GENE = mod2gene, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
height_modGSEA <- height_modGSEA@result
# alzheimer_modGSEA <- GSEA(alzheimer_snpgenelist, TERM2GENE = mod2gene, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
# alzheimer_modGSEA <- alzheimer_modGSEA@result

diab_submodGSEA001 <- GSEA(diab_snpgenelist, TERM2GENE = valuable_subMod001, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
diab_submodGSEA001 <- diab_submodGSEA001@result
height_submodGSEA001 <- GSEA(height_snpgenelist, TERM2GENE = valuable_subMod001, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
height_submodGSEA001 <- height_submodGSEA001@result
# alzheimer_submodGSEA001 <- GSEA(alzheimer_snpgenelist, TERM2GENE = valuable_subMod001, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
# alzheimer_submodGSEA001 <- alzheimer_submodGSEA001@result

diab_submodGSEA01 <- GSEA(diab_snpgenelist, TERM2GENE = valuable_subMod01, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
diab_submodGSEA01 <- diab_submodGSEA01@result
height_submodGSEA01 <- GSEA(height_snpgenelist, TERM2GENE = valuable_subMod01, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
height_submodGSEA01 <- height_submodGSEA01@result
# alzheimer_submodGSEA01 <- GSEA(alzheimer_snpgenelist, TERM2GENE = valuable_subMod01, pvalueCutoff = 1, minGSSize = 0, maxGSSize = 10000, nPerm = 20000, exponent = 1)
# alzheimer_submodGSEA01 <- alzheimer_submodGSEA01@result

## write out the GSEA results
tmp <- data.frame(diab_modGSEA[, c("ID", "pvalue")], height_modGSEA[, c("ID", "pvalue")])
colnames(tmp) <- c("Diab_ID", "Diab_p", "height_ID", "height_p")
write.table(tmp, file = "./caudate/mod_GSEA.txt", quote = F, sep = "\t", row.names = F)
tmp <- data.frame(diab_submodGSEA01[, c("ID", "pvalue")], height_submodGSEA01[, c("ID", "pvalue")])
colnames(tmp) <- c("Diab_ID", "Diab_p", "height_ID", "height_p")
write.table(tmp, file = "./caudate/submod01_GSEA.txt", quote = F, sep = "\t", row.names = F)


# #--------------------------------------------------------------------------------------------------------------
# # IV. mod GO
# #--------------------------------------------------------------------------------------------------------------
# GOmodgene <- list(Turquoise = modgene$turquoise, Down = downgene, Up = upgene, SUB = intersect(modgene$turquoise, downgene))
# GOsubmodgene <- split(valuable_subMod01$Gene, valuable_subMod01$subMod)
# GOsubmodgene <- GOsubmodgene[c("MEturquoise_Cpu.D2_0.01", "MEturquoise_Cpu.D1_0.01", "MEturquoise_Striatum.Young.Adulthood_0.01")]
# 
# ## GO function    
# modGO <- function(genelist, keyType) {
#     if (keyType == "ENSEMBL") {
#         genelist <- sub("\\.[0-9]+", "", genelist)
#         readable <-  TRUE
#     } else if (keyType != "SYMBOL") {
#         message("Can't recognize this keyType!")
#     } else {readable <-  FALSE}
#     ego <- enrichGO(gene = genelist,
#                     OrgDb = org.Hs.eg.db,
#                     keyType = keyType,
#                     ont = "BP",
#                     pvalueCutoff = 0.01,
#                     qvalueCutoff = 0.05,
#                     readable = readable)
#     data.frame(ego)
# }
# 
# ## KEGG function 
# modKEGG <- function(genelist, keyType) {
#     if (keyType == "ENSEMBL") {
#         genelist <- sub("\\.[0-9]+", "", genelist)
#         genelist <- bitr(genelist, fromType = "ENSEMBL", toType = "UNIPROT", OrgDb = "org.Hs.eg.db")
#     } else if (keyType == "SYMBOL") {
#         genelist <- bitr(genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = "org.Hs.eg.db")
#     } else {
#         message("Can't recognize this keyType!")
#     }
#     kk <- enrichKEGG(gene = genelist$UNIPROT,
#                      organism = "hsa",
#                      keyType = "uniprot",
#                      pvalueCutoff = 0.05)
#     data.frame(kk)
# }
# 
# GOres <- lapply(X = GOmodgene, FUN = modGO, "ENSEMBL")
# KEGGres <- lapply(X = GOmodgene, FUN = modKEGG, "ENSEMBL")
# GOsubres <- lapply(X = GOsubmodgene, FUN = modGO, "SYMBOL")
# KEGGsubres <- lapply(X = GOsubmodgene, FUN = modKEGG, "SYMBOL")
# 
# ## convert list to dataframe
# list2data.frame <- function(x) {
#     listmodname <- names(x)
#     res_df <- data.frame(stringsAsFactors = F)
#     for (i in 1:length(x)) {
#         tmp <- x[[i]]
#         tmp <- data.frame(Mod = rep(listmodname[i], dim(tmp)[1]), tmp, stringsAsFactors = F)
#         res_df <- rbind(res_df, tmp)
#     }
#     res_df
# }
# 
# 
# GOdf <- list2data.frame(GOres)
# KEGGdf <- list2data.frame(KEGGres)
# GOsubdf <- list2data.frame(GOsubres)
# KEGGsubdf <- list2data.frame(KEGGsubres)
# write.table(GOdf, file = "./caudate/GO.tab", quote = F, sep = "\t", row.names = F)
# write.table(KEGGdf, file = "./caudate/KEGG.tab", quote = F, sep = "\t", row.names = F)
# 
# #save(GOdf, KEGGdf, mod_summary, mod2gene, modgene, DEgene, bwnet, file = "../caudate_covar_WGCNA.RData")
