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
rm(list=ls())
options(stringsAsFactors = FALSE)

snp_enrich_plot <- function(permutation_data, mod, color, gwas_type) {
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
    hist((permutation_ES-meanx)/sdx, breaks = 100, xlim = c(-7, 7), 
         col = color, main = gwas_type, xlab = paste0("Z-score of ES"))
    axis(1, at = c(-7:7))
    abline(v = (true_ES-meanx)/sdx, col = "red")
}

# ----------------------------------------------------------------------------------------------------
# I. caudate
# ----------------------------------------------------------------------------------------------------

inputdir <- "/data/MyProgram/Final_diabrain/7.scripts/7_GWAS_REGION/script/caudate/"
setwd(inputdir)
outdir <- "/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/"
pdf(paste0(outdir, "GWAS_plot_all_module.pdf"), width = 12, height = 7)

# split the canvas and prepare a out margin for title
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))

# plot the two hist
diab <- read.table("./diab/edb/results.txt", sep = "\t")
height <- read.table("./height/edb/results.txt", sep = "\t")

mname <- tolower(sub("ME", "", diab$V1))
for (module in mname) {
    if (module %in% c("upgene", "downgene", "totalgene")) {mycolor <- "grey"} else {mycolor <- module}
    snp_enrich_plot(diab, mod = module, color = mycolor, gwas_type = "T2D")
    snp_enrich_plot(height, mod = module, color =mycolor, gwas_type = "Height")
    # title
    mtext(paste("SNP enrichment in", module, "module"), side = 3, line = 0, outer = T)
}

dev.off()


# ----------------------------------------------------------------------------------------------------
# I. sub_mod
# ----------------------------------------------------------------------------------------------------
diab_sub <- read.table("./diab_sub.GseaPreranked.1552799465391/edb/results.txt", sep = "\t")
height_sub <- read.table("./height_sub.GseaPreranked.1552799533970/edb/results.txt", sep = "\t")

pdf(paste0(outdir, "GWAS_plot_sub_module.pdf"), width = 12, height = 7)
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
