# ====================================================================================================
# Author:   Aaron
# Function: plot the permutaions of GWAS gene GSEA results
# Version:  1.0
# Date:     Mar 18, 2019
# ====================================================================================================
#library(clusterProfiler)
#library(DOSE)
#library(org.Hs.eg.db)
#library(WGCNA)
library(extrafont)
loadfonts()
rm(list=ls())
options(stringsAsFactors = FALSE)

snp_enrich_plot <- function(permutation_data, mod, color, gwas_type, main = NULL) {
    mod_name <- tolower(sub("ME", "", permutation_data$V1))
    diab_esmatrix <- permutation_data$V13
    names(diab_esmatrix) <- mod_name
    diab_esmatrix <- strsplit(diab_esmatrix, split = " ")
    diab_esmatrix <- sapply(diab_esmatrix, function(x) {as.numeric(x)})
    
    permutation_ES <- diab_esmatrix[, mod]
    true_ES <- permutation_data$V3[mod_name == mod]
    
    meanx <- mean(permutation_ES)
    sdx <- sd(permutation_ES)
    p_value <- signif(sum(permutation_ES > true_ES)/20000, digits = 4)
    zdata <- (permutation_ES-meanx)/sdx
    maxn <- ceiling(max(abs(zdata)))
    if (is.null(main)) {main <- gwas_type}
    hist(zdata, breaks = seq(-maxn, maxn, 0.2), xlim = c(-7, 7), 
         col = color, main = main, xlab = paste0("Z-score of ES"), 
         cex.lab = 1.8, cex.main = 1.8, cex.axis = 1.1, lwd = 1.5)
    axis(1, at = c(-7:7), cex.axis = 1.1, lwd = 1.5)
    abline(v = (true_ES-meanx)/sdx, col = "red")
}

# ----------------------------------------------------------------------------------------------------
# I. caudate
# ----------------------------------------------------------------------------------------------------

inputdir <- "/data/MyProgram/Final_diabrain/7.interdata/7_GWAS_REGION/script/caudate/"
setwd(inputdir)
outdir <- "/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/"
pdf(paste0(outdir, "GWAS_plot_all_module.pdf"), width = 12, height = 5)

# split the canvas and prepare a out margin for title
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))

# plot the two hist
diab <- read.table("./diab/edb/results.txt", sep = "\t")
height <- read.table("./height/edb/results.txt", sep = "\t")

mname <- tolower(sub("ME", "", diab$V1))
for (module in mname) {
    # if (module %in% c("upgene", "downgene", "totalgene")) {mycolor <- "grey"} else {mycolor <- module}
    if (module %in% c("upgene", "downgene", "totalgene")) {next} else {mycolor <- module}
    snp_enrich_plot(diab, mod = module, color = mycolor, gwas_type = "T2D")
    snp_enrich_plot(height, mod = module, color =mycolor, gwas_type = "Height")
    # title
    mtext(paste("SNP enrichment in", module, "module"), side = 3, line = 0, outer = T)
}

dev.off()


# ----------------------------------------------------------------------------------------------------
# II. sub_mod
# ----------------------------------------------------------------------------------------------------
diab_sub <- read.table("./diab_sub.GseaPreranked.1552799465391/edb/results.txt", sep = "\t")
height_sub <- read.table("./height_sub.GseaPreranked.1552799533970/edb/results.txt", sep = "\t")

pdf(paste0(outdir, "GWAS_plot_sub_module.pdf"), family = "Arial", width = 12, height = 5)
# cpu.D2
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
snp_enrich_plot(diab_sub, mod = "turquoise_cpu.d2_0.01", color = "turquoise", gwas_type = "T2D")
snp_enrich_plot(height_sub, mod = "turquoise_cpu.d2_0.01", color = "turquoise", gwas_type = "Height")
mtext("SNP enrichment in D2-MSN sub-module", side = 3, line = 0, outer = T)

# cpu.D1
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
snp_enrich_plot(diab_sub, mod = "turquoise_cpu.d1_0.01", color = "turquoise", gwas_type = "T2D")
snp_enrich_plot(height_sub, mod = "turquoise_cpu.d1_0.01", color = "turquoise", gwas_type = "Height")
mtext("SNP enrichment in D1-MSN sub-module", side = 3, line = 0, outer = T)

# striatum
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
snp_enrich_plot(diab_sub, mod = "turquoise_striatum.young.adulthood_0.01", color = "turquoise", gwas_type = "T2D")
snp_enrich_plot(height_sub, mod = "turquoise_striatum.young.adulthood_0.01", color = "turquoise", gwas_type = "Height")
mtext("SNP enrichment in striatum sub-module", side = 3, line = 0, outer = T)

dev.off()
# embed("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/GWAS_plot_sub_module.pdf")

# ----------------------------------------------------------------------------------------------------
# III. interest
# ----------------------------------------------------------------------------------------------------
pdf(paste0(outdir, "GWAS_plot_interest_mod_new.pdf"), width = 12, height = 6)
par(mfrow = c(2, 4), ps=12)
for (module in c("turquoise", "red", "magenta", "green", "tan")) {
    par(mar=c(5,5,4,0.5))
    snp_enrich_plot(diab, mod = module, color = module, gwas_type = "T2D", main = module)
}
snp_enrich_plot(diab_sub, mod = "turquoise_cpu.d2_0.01", color = "turquoise", gwas_type = "T2D", main = paste("MSN-D2"))
snp_enrich_plot(diab_sub, mod = "turquoise_cpu.d1_0.01", color = "turquoise", gwas_type = "T2D", main = paste("MSN-D1"))
snp_enrich_plot(diab_sub, mod = "turquoise_striatum.young.adulthood_0.01", color = "turquoise", gwas_type = "T2D",  main = paste("striatum"))
dev.off()

# plot the pvalue of interest mod
# p.df <- data.frame(Mod = c("turquoise", "red", "magenta", "green", "tan", "MSN-D2"), 
#                    T2D_p = c(0, 0.00815, 0.0091, 0.0095, 0.01205, 0.0144), 
#                    Height_p = c(0.0497, 0.2414, 0.1245, 0.04135, 0.24836242, 0.45245))
p.df <- read.table("./Table S8.csv", header = T, sep = ",")
p.df <- p.df[c(1:5, 15:17), c("Module", "T2D.FDR", "Height.FDR")]
for_plot <- data.frame(logp = log10(c(p.df$T2D.FDR, p.df$Height.FDR)), Type = rep(c("T2D", "Height"), each = dim(p.df)[1]), module = rep(p.df$Module, 2))
for_plot$module  <- factor(for_plot$module, levels = c("turquoise", "red", "magenta", "green", "tan", "D2-MSN", "D1-MSN", "striatum"))
# make the infinity to a relative large value
for_plot[1,1] <- -5

pdf(paste0(outdir, "GWAS_pvalue.pdf"), family = "Arial", width = 10, height = 6)
library(ggplot2)
ggplot(for_plot, aes(x=module,  y = logp, group = Type, color = Type, fill = Type)) + 
    geom_bar(stat = "identity", width = 0.5, position = "dodge") + 
    geom_hline(aes(yintercept=log10(0.05)), colour = "red") +
    scale_fill_manual(values = c("grey", "red")) +
    scale_color_manual(values = c("grey", "red")) +
    # coord_flip() + 
    labs(x="Mod",y="P-value")+
    scale_y_continuous(breaks = c(-5:0), labels = c("0", "0.0001", "0.001", "0.01", "0.1", "1"), expand = c(0, 0)) + 
    scale_x_discrete(position = "top") +
    theme(axis.title.x=element_text(size=15,face="bold"),
          axis.title.y=element_text(size=15,face="bold"),
          axis.text=element_text(size=15),
          axis.text.x = element_text(face = "bold"),
          legend.title=element_text(size=15,face="bold"),
          legend.text=element_text(size=15,face="bold"),
          panel.grid.major  = element_blank(),
          panel.background = element_blank(), 
          axis.line=element_line(size = .6, colour="black"), 
          strip.text.x=element_text(size=15,face="bold"))
    
# abline(h = log(0.05))
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/GWAS_pvalue.pdf")
