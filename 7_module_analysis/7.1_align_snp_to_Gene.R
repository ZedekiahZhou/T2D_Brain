
# ====================================================================================================
# Author:   Aaron
# Function: map snp to gene
# Version:  1.0
# Date:     Jan 9, 2019
# ====================================================================================================
# library(biomaRt)
rm(list=ls())
options(stringsAsFactors = FALSE)
setwd("/data/MyProgram/Final_diabrain/7.scripts/7_GWAS_REGION/script/")
auto_gene <- read.table("1_Gene_anno_autosome.tab", stringsAsFactors = F, header = T)
chrlen <- read.table("1_GRCh37_chr_len.tab", stringsAsFactors = F, header = T)
diab_snp <- read.table("0_Xue_DIABETES_snp_summary_short.txt", stringsAsFactors = F, header = T)
height_snp <- read.table("0_height_standart_snp.txt", stringsAsFactors = F, header = T)
alzheimer_snp <- read.table("0_Jansen_Alzheimer_snp_summary_short.txt", header = T)

# modify some attributes to fit the functions
height_snp$CHR <- as.integer(height_snp$CHR)
# names(alzheimer_snp) <- c("CHR", "BP", "SNP", "P")

diab_snp <- diab_snp[order(diab_snp$CHR, diab_snp$BP), ]
height_snp <- height_snp[order(height_snp$CHR, height_snp$BP), ]
alzheimer_snp <- alzheimer_snp[order(alzheimer_snp$CHR, alzheimer_snp$BP), ]
auto_gene <- auto_gene[order(auto_gene$chr, auto_gene$start, auto_gene$end), ]

# --------------------------------------------------------------------------
# the function to judge if the snp locus whin one gene region
# --------------------------------------------------------------------------
snpingene <- function(snp_chr, snp_pos, gene_chr, gene_up, gene_down) {
    if (snp_chr == gene_chr) {
        if (snp_pos < gene_up) {
            return("less")
        } else if (snp_pos > gene_down) {
            return("greater")
        } else {
            return("yes")
        }
    } else if (snp_chr < gene_chr) {
        return("less")
    } else {
        return("greater")
    }
}

# ---------------------------------------------------------------------------------------------------------------
# the function to map snp to gene
# ---------------------------------------------------------------------------------------------------------------
bi_mapSNPtoGENE <- function(snp, gene) {
    gene_stas <- data.frame(SNP = rep("Nosnp", dim(gene)[1]), pvalue = rep(1, dim(gene)[1]), stringsAsFactors = F)
    
    # Method 1: Bisection query
    n1 <- 0
    for (i in 1:dim(gene)[1]) {
        lpoint <- 1; rpoint <- dim(snp)[1]
        genel <- gene$start[i] - 20000
        gener <- gene$end[i] + 10000
        genec <- gene$chr[i]
        if ((snpingene(snp$CHR[lpoint], snp$BP[lpoint], genec, genel, gener)=="greater") | 
            (snpingene(snp$CHR[rpoint], snp$BP[rpoint], genec, genel, gener)=="less")) next
        
        while (lpoint <= rpoint) {
            mi <- ceiling((lpoint+rpoint)/2)
            idx <- snpingene(snp$CHR[mi], snp$BP[mi], genec, genel, gener)
            if( idx == "less") {
                lpoint <- mi + 1
                n1 <- n1+1
            } else if (idx == "greater") {
                rpoint <- mi - 1
                n1 <- n1+1
            } else if (idx == "yes") {
                x <- mi
                while ((snp$CHR[x] == genec) & (snp$BP[x] >= genel)) {
                    if (snp$P[x] < gene_stas$pvalue[i]) {
                        gene_stas$SNP[i] <- snp$SNP[x]
                        gene_stas$pvalue[i] <- snp$P[x]
                    }
                    x <- x-1
                    n1 <- n1 +1
                    if (x < 1) break
                }
                x <- mi+1
                while ((snp$CHR[x] == genec) & (snp$BP[x] <= gener)) {
                    if (snp$P[x] < gene_stas$pvalue[i]) {
                        gene_stas$SNP[i] <- snp$SNP[x]
                        gene_stas$pvalue[i] <- snp$P[x]
                    }
                    x <- x+1
                    n1 <- n1 + 1
                    if (x > dim(snp)[1]) break
                }
                lpoint <- rpoint + 1
            }
        }
    }
    gene_snp <- cbind(gene, gene_stas)
    # gene_snp <- subset(gene, SNP!= "Nosnp")
    
    gene_snp <- cbind(gene_snp, minus_log10P = -log10(gene_snp$pvalue))
    gene_snp <- gene_snp[order(gene_snp$minus_log10P, decreasing = T), ]
    return(gene_snp)
}

system.time(diab_snp_to_gene <- bi_mapSNPtoGENE(diab_snp, auto_gene))
system.time(height_snp_to_gene <- bi_mapSNPtoGENE(height_snp, auto_gene))
system.time(alzheimer_snp_to_gene <- bi_mapSNPtoGENE(alzheimer_snp, auto_gene))

write.table(diab_snp_to_gene, file = "2_DIABETES_SNP_GENEs_for_GSEA.tab", quote = F, sep = "\t", row.names = F)
write.table(height_snp_to_gene, file = "2_HEIGHT_SNP_GENEs_for_GSEA.tab", quote = F, sep = "\t", row.names = F)
write.table(alzheimer_snp_to_gene, file = "2_ALZHEIMER_SNP_GENES_for_GSEA.tab", quote = F, sep = "\t", row.names = F)


