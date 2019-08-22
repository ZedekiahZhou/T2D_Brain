# ====================================================================================================
# Author:   Aaron
# Function: Pathway analysis
#           Prepare genes for TF analysis
# Version:  1.0
# Date:     Mar 19, 2019
# ====================================================================================================

library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(WGCNA)
library(extrafont)
loadfonts()

#--------------------------------------------------------------------------------------------------------------
# I. mod GO
#--------------------------------------------------------------------------------------------------------------
rm(list=ls())
setwd("/data/MyProgram/Final_diabrain/4.plots/GO/")
options(stringsAsFactors = FALSE)

# load mod information
tmp <- load("/data/MyProgram/Final_diabrain/7.interdata/6_WGCNA/Brain_Caudate_covar_WGCNAres_v1.0.RData")
tmp <- load("/data/MyProgram/Final_diabrain/7.interdata/7_GWAS_REGION/script/caudate/region_v1.0.RData")

## GO & KEGG function 
modPathway <- function(genelist, keyType, outdir, genelist.name) {
    if (keyType == "ENSEMBL") {
        genelist <- sub("\\.[0-9]+", "", genelist)
        prolist <- bitr(genelist, fromType = "ENSEMBL", toType = "UNIPROT", OrgDb = "org.Hs.eg.db")
        readable <-  TRUE
    } else if (keyType == "SYMBOL") {
        prolist <- bitr(genelist, fromType = "SYMBOL", toType = "UNIPROT", OrgDb = "org.Hs.eg.db")
        readable <-  FALSE
    } else {
        message("Can't recognize this keyType!")
    }
    ego <- enrichGO(gene = genelist,
                    OrgDb = org.Hs.eg.db,
                    keyType = keyType,
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    readable = readable)
    # ego_simp <- simplify(ego)
    egodf <- data.frame(ego)
    kk <- enrichKEGG(gene = prolist$UNIPROT,
                     organism = "hsa",
                     keyType = "uniprot",
                     pvalueCutoff = 0.05)
    kkdf <- data.frame(kk)
    save(ego, kk, file = paste0(outdir, genelist.name, ".RData"))
    rbind(egodf, kkdf)
}

modgene <- split(names(datExpr), f = as.factor(bwnetColors))
# turquoisedf <- modPathway(genelist = modgene$turquoise, keyType = "ENSEMBL", outdir = "./caudate_mod_ redunt/", genelist.name = "turquoise")

dir.create("./caudate_mod_ redunt/")
allmoddf <- list()
for (i in 1:length(modgene)) {
    allmoddf[[i]] <- modPathway(modgene[[i]], keyType = "ENSEMBL", outdir = "./caudate_mod_ redunt/", genelist.name = names(modgene)[i])
    write.table(allmoddf[[i]], paste0("./caudate_mod_ redunt/", names(modgene)[i], ".tab"), quote = F, sep = "\t", row.names = F)
}

# module go analysis
GOsubmodgene <- split(valuable_subMod01$Gene, valuable_subMod01$subMod)
D1df <- modPathway(genelist = GOsubmodgene$MEturquoise_Cpu.D1_0.01, keyType = "SYMBOL", outdir = "./", genelist.name = "Caudate_D1MSN")
D2df <- modPathway(genelist = GOsubmodgene$MEturquoise_Cpu.D2_0.01, keyType = "SYMBOL", outdir = "./", genelist.name = "Caudate_D2MSN")
Striatumdf <- modPathway(genelist = GOsubmodgene$MEturquoise_Striatum.Young.Adulthood_0.01, keyType = "SYMBOL", outdir = "./", genelist.name = "Caudate_striatum")

# write out the GO results
# write.table(turquoisedf, "./Caudate_turquoise.tab", quote = F, sep = "\t", row.names = F)
write.table(D1df, "./Caudate_D1MSN.tab", quote = F, sep = "\t", row.names = F)
write.table(D2df, "./Caudate_D2MSN.tab", quote = F, sep = "\t", row.names = F)
write.table(Striatumdf, "./Caudate_striatum.tab", quote = F, sep = "\t", row.names = F)

# write out the genes
write(sub("\\.[0-9]+", "", modgene$turquoise), "../TF/Caudate_turguoise_gene.list")
write(GOsubmodgene$MEturquoise_Cpu.D1_0.01, "../TF/Caudate_D1MSN_gene.list")
write(GOsubmodgene$MEturquoise_Cpu.D2_0.01, "../TF/Caudate_D2MSN_gene.list")
write(GOsubmodgene$MEturquoise_Striatum.Young.Adulthood_0.01, "../TF/Caudate_Striatum_gene.list")

#--------------------------------------------------------------------------------------------------------------
# II. DEGs GO
#--------------------------------------------------------------------------------------------------------------
Caudate_up <- read.table("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/result05/Up_Gene/Brain_Caudate.up")[[1]]
Caudate_down <- read.table("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/result05/Down_Gene/Brain_Caudate.down")[[1]]
Caudate_up_df <- modPathway(genelist = Caudate_up, keyType = "ENSEMBL", outdir = "./", genelist.name = "Caudate_up")
Caudate_down_df <- modPathway(genelist = Caudate_down, keyType = "ENSEMBL", outdir = "./", genelist.name = "Caudate_down")

Hippocampus_up <- read.table("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/result05/Up_Gene/Brain_Hippocampus.up")[[1]]
Hippocampus_down <- read.table("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/result05/Down_Gene/Brain_Hippocampus.down")[[1]]
Hippocampus_up_df <- modPathway(genelist = Hippocampus_up, keyType = "ENSEMBL", outdir = "./", genelist.name = "Hippocampus_up")
Hippocampus_down_df <- modPathway(genelist = Hippocampus_down, keyType = "ENSEMBL", outdir = "./", genelist.name = "Hippocampus_down")

# write out go res
write.table(Caudate_up_df, "./Caudate_up.tab", quote = F, sep = "\t", row.names = F)
write.table(Caudate_down_df, "./Caudate_down.tab", quote = F, sep = "\t", row.names = F)
write.table(Hippocampus_up_df, "./Hippocampus_up.tab", quote = F, sep = "\t", row.names = F)
write.table(Hippocampus_down_df, "./Hippocampus_down.tab", quote = F, sep = "\t", row.names = F)

# write out the genes
write(sub("\\.[0-9]+", "", Caudate_up), "../TF/Caudate_up_gene.list")
write(sub("\\.[0-9]+", "", Caudate_down), "../TF/Caudate_down_gene.list")
write(sub("\\.[0-9]+", "", Hippocampus_up), "../TF/Hippocampus_up_gene.list")
write(sub("\\.[0-9]+", "", Hippocampus_down), "../TF/Hippocampus_down_gene.list")

#--------------------------------------------------------------------------------------------------------------
# III. GO plot for turquoise
#--------------------------------------------------------------------------------------------------------------
turquoise.anno <- read.table("/data/MyProgram/Final_diabrain/4.plots/GO/caudate_mod_ redunt/turquoise.tab", sep = "\t", header = T)
turquoise.anno$Class <- gsub("[[:digit:]]+", "", turquoise.anno$ID)
turquoise.anno$Name <- paste0(turquoise.anno$Description, "(", turquoise.anno$Count, ")")
turquoise.anno$minus.log10p <- -log10(turquoise.anno$p.adjust)
tmp <- split(turquoise.anno, f = turquoise.anno$Class)

# construct predata for ggplot2
go.pre <- tmp$"GO:"[1:20, c("Name", "minus.log10p")]
kegg.pre <- tmp$hsa[1:20, c("Name", "minus.log10p")]
# specify the order
go.pre$Name <- factor(go.pre$Name, levels = rev(go.pre$Name))
kegg.pre$Name <- factor(kegg.pre$Name, levels = rev(kegg.pre$Name))

library(ggplot2)
library(grid)
library(RColorBrewer)
mycol <- brewer.pal(11, "Spectral")[10]
#mycol <- "turquoise"
p.go <- ggplot(go.pre, aes(x=Name, y = minus.log10p)) + 
    geom_bar(width = 0.5, stat = "identity", fill = mycol, colour = mycol) + 
    coord_flip() + 
    labs(x="",y="-log10(adjust P-value)")+
    scale_y_continuous(expand = c(0, 0)) +   # remove the gap between geom and axis
    theme(axis.title=element_text(size=8,face="bold"),
          axis.text=element_text(size=8),
          axis.text.x.top = element_text(face = "bold"),
          legend.key.height = unit(8, "points"), 
          legend.key.width = unit(10, "points"),
          legend.position = "top",
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(face="bold"),
          panel.grid.major  = element_blank(),     # remove the grid
          panel.background = element_blank(),      # remove the grey backgroud
          axis.line=element_line(size = .3, colour="black"), 
          strip.text.x=element_text(size=8,face="bold"))
p.kegg <- ggplot(kegg.pre, aes(x=Name, y = minus.log10p)) + 
    geom_bar(width = 0.5, stat = "identity", fill = mycol, colour = mycol) + 
    coord_flip() + 
    labs(x="",y="-log10(adjust P-value)")+
    scale_y_continuous(expand = c(0, 0)) +   # remove the gap between geom and axis
    theme(axis.title=element_text(size=8,face="bold"),
          axis.text=element_text(size=8),
          axis.text.x.top = element_text(face = "bold"),
          legend.key.height = unit(8, "points"), 
          legend.key.width = unit(10, "points"),
          legend.position = "top",
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(face="bold"),
          panel.grid.major  = element_blank(),     # remove the grid
          panel.background = element_blank(),      # remove the grey backgroud
          axis.line=element_line(size = .3, colour="black"), 
          strip.text.x=element_text(size=8,face="bold"))

pdf(file = "./turquoise.pathway.plot.blue.pdf", family = "Arial", width = 8, height = 4)
pushViewport(viewport(layout = grid.layout(1,2)))     # layout
print(p.go, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p.kegg, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
