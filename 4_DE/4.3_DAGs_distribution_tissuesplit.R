# ====================================================================================================
# Author:   Aaron
# Function: DEGs distribution for two tissue
# Version:  1.0
# Date:     May 7, 2019
# ====================================================================================================
library(ggplot2)
library(ggpubr)
library(scales)
library(extrafont)
loadfonts()
library(RColorBrewer)
rm(list = ls())
options(stringsAsFactors = FALSE)

inputfiles <- dir("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/result05/res/", full.names = T)
tnames <- sub(".res", "", sub("[[:print:]]+Brain_", "", inputfiles))
geneanno <- read.table("/data/MyProgram/Final_diabrain/1.clean/gen_anno/Total_Gene_Info2.list", header = T)
mycolor <- brewer.pal(11, "Spectral")[c(2,10, 11)]


# ----------------------------------------------------------------------------------------------------
# I. pooled the DEGs from all the 13 tissue
# ----------------------------------------------------------------------------------------------------
allDAGs <- data.frame()
for (i in 1:length(inputfiles)) {
    tmp <- read.table(inputfiles[i], header = T)
    if (nrow(tmp)==0) next
    tmp <- cbind(tmp, REGION = rep(tnames[i], nrow(tmp)), geneanno[match(tmp$Gene, geneanno$ID), ])
    allDAGs <- rbind(allDAGs, tmp)
}
write.table(allDAGs, "/data/MyProgram/Final_diabrain/4.plots/DEGs/DAGs_info.txt", sep = "\t", quote = F, row.names = F)


# get the annotations of all genes
total <- cbind(allDAGs, Regulation = ifelse(allDAGs$log2FoldChange > 0, "Up", "Down"))[, c("Symbol", "REGION", "GeneType", "Regulation", "CHR", "Start", "End")]
caudate <- subset(total, REGION == "Caudate")
# caudate_up <- subset(total, REGION == "Caudate" & Regulation == "Up")
# caudate_down <- subset(total, REGION == "Caudate" & Regulation == "Down")
hippocampus <- subset(total, REGION == "Hippocampus")
# hippocampus_up <- subset(total, REGION == "Hippocampus" & Regulation == "Up")
# hippocampus_down <- subset(total, REGION == "Hippocampus" & Regulation == "Down")

# ----------------------------------------------------------------------------------------------------
# II. Gene type
# ----------------------------------------------------------------------------------------------------
# prepare for bar plot
pre_data <- rbind(caudate, hippocampus)

# summary for the genetype counts
cau_gene_type <- as.data.frame(table(pre_data$GeneType[pre_data$REGION == "Caudate"]))
hippo_gene_type <- as.data.frame(table(pre_data$GeneType[pre_data$REGION == "Hippocampus"]))
colnames(cau_gene_type) <- c("GeneType", "Freq.Caudate")
colnames(hippo_gene_type) <- c("GeneType", "Freq.Hippocampus")
cau_gene_type$Percent.Caudate <- percent(cau_gene_type$Freq.Caudate/sum(cau_gene_type$Freq.Caudate))
hippo_gene_type$Percent.Hippocampus <- percent(hippo_gene_type$Freq.Hippocampus/sum(hippo_gene_type$Freq.Hippocampus))
total_type <- merge(cau_gene_type, hippo_gene_type, by = "GeneType")
write.table(total_type, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Gene_type_distribution_tissuesplit.txt", sep = "\t", quote = F, row.names = F)

genetype_description <- data.frame(Genetype = c("protein_coding", "lincRNA", "antisense", "pseudogene"), 
                                   Description = c("Protein coding", "LincRNA", "Antisense RNA", "Pseudogene"))
idx <- match(pre_data$GeneType, genetype_description$Genetype)
pre_data$GeneType <- ifelse(is.na(idx), "Others", genetype_description$Description[idx])
pre_data$GeneType <- as.factor(pre_data$GeneType)

# bar plot
pdf("/data/MyProgram/Final_diabrain/4.plots/DEGs/Gene_type_distribution_tissuesplit.pdf", family = "Arial", height = 2.5, width = 4)
ggplot(pre_data, aes(x=GeneType, y = ..prop.., group = REGION, color = REGION, fill = REGION)) + 
    geom_bar(width = 0.7, position = position_dodge2(reverse = TRUE)) +      # plot groups side by side
    coord_flip() +                 # reverse x and y axis
    scale_y_continuous(labels = percent, expand = c(0, 0)) +    # remove the gab between geom and x axis
    scale_fill_manual(values = mycolor[1:2]) + 
    scale_color_manual(values = mycolor[1:2]) +
    labs(x="Gene Type",y="Proportion")+
    theme(axis.title=element_text(size=8,face="bold"),
          axis.text=element_text(size=8),
          axis.text.y = element_text(face = "bold"),
          legend.title=element_text(size=8,face="bold"),
          legend.key.height = unit(8, "points"), 
          legend.key.width = unit(10, "points"),
          legend.position = "top",
          # legend.title=element_blank(),
          legend.text=element_text(face="bold"),
          panel.grid.major  = element_blank(),            # remove backgroud 
          panel.background = element_blank(),             # remove backgroud 
          axis.line=element_line(size = .3, colour="black"),  # add axis
          strip.text.x=element_text(size=8,face="bold"))
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/DEGs/Gene_type_distribution_tissuesplit.pdf")



# ----------------------------------------------------------------------------------------------------
# III. chromosome distribution
# ----------------------------------------------------------------------------------------------------
pre_data$CHR <- factor(pre_data$CHR, levels = c(as.character(1:22), "X", "Y", "MT"))

# summary for the chromosome counts
cau_gene_chr <- as.data.frame(table(pre_data$CHR[pre_data$REGION == "Caudate"]))
hippo_gene_chr <- as.data.frame(table(pre_data$CHR[pre_data$REGION == "Hippocampus"]))
colnames(cau_gene_chr) <- c("Chr", "Freq.Caudate")
colnames(hippo_gene_chr) <- c("Chr", "Freq.Hippocampus")
cau_gene_chr$Percent.Caudate <- percent(cau_gene_chr$Freq.Caudate/sum(cau_gene_chr$Freq.Caudate))
hippo_gene_chr$Percent.Down <- percent(hippo_gene_chr$Freq.Hippocampus/sum(hippo_gene_chr$Freq.Hippocampus))
total_chr <- merge(cau_gene_chr, hippo_gene_chr, by = "Chr")
total_chr <- total_chr[order(total_chr$Chr), ]
write.table(total_chr, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distribution_tissuesplit.txt", sep = "\t", quote = F, row.names = F)

pdf("/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distribution_tissuesplit.pdf", family = "Arial", height = 2.5, width = 8)
ggplot(pre_data, aes(x=CHR, y = ..prop.., group = REGION, color = REGION, fill = REGION)) + 
    geom_bar(width = 0.5, position = "dodge") + 
    # coord_flip() + 
    labs(x="Chromosome",y="Proportion")+
    scale_y_continuous(labels = percent, expand = c(0, 0)) + 
    scale_fill_manual(values = mycolor[1:2]) +       # specify the color
    scale_color_manual(values = mycolor[1:2]) +
    theme(axis.title=element_text(size=8,face="bold"),
          axis.text=element_text(size=8),
          axis.text.x = element_text(face = "bold"),
          legend.key.height = unit(8, "points"), 
          legend.key.width = unit(10, "points"),
          legend.position = "top",
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(face="bold"),
          panel.grid.major  = element_blank(),
          panel.background = element_blank(), 
          axis.line=element_line(size = .3, colour="black"), 
          strip.text.x=element_text(size=8,face="bold"))
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distribution_tissuesplit.pdf")

# ----------------------------------------------------------------------------------------------------
# IV. chromosome distance between region
# ----------------------------------------------------------------------------------------------------
# pre_data <- pre_data[order(pre_data$CHR), ]
chr_split <- split(pre_data, f = pre_data$CHR)
chr_dist = as.integer()
Gene_Pair = as.character()
Pair_Type = as.character()
Region_Type <- as.character() -> Regulation_Type
l <- 0

for (i in 1:length(chr_split)) {
    if (dim(chr_split[[i]])[1] == 1) next
    for (j in 1:(dim(chr_split[[i]])[1]-1) ) {
        for (k in (j+1):dim(chr_split[[i]])[1] ) {
            dis1 <- chr_split[[i]]$Start[j] - chr_split[[i]]$End[k]
            dis2 <- chr_split[[i]]$Start[k] - chr_split[[i]]$End[j]
            l <- l+1
            chr_dist[l] <- max(dis1, dis2, 0)
            Gene_Pair[l] <- paste(chr_split[[i]]$Symbol[j], chr_split[[i]]$Symbol[k], sep = "-")
            if (chr_split[[i]]$REGION[j] != chr_split[[i]]$REGION[k]) {
                Region_Type[l] <- "CHs"
            } else if (chr_split[[i]]$REGION[j] == "Caudate") {
                Region_Type[l] <- "CCs"
            } else {
                Region_Type[l] <- "HHs"
            }
            # Regulation_Type[l] <- paste(chr_split[[i]]$Regulation[j], chr_split[[i]]$Regulation[k], sep = "-")
        }
    }
}


# Regulation_Type[Regulation_Type == "Down-Up"] <- "Up-Down"
chrdis_df <- data.frame(chr_dist, Gene_Pair, Region_Type)
chrdis_df$Region_Type <- factor(chrdis_df$Region_Type, levels = c("CCs", "HHs", "CHs"))

pdf("/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distance_tissuesplit.pdf", family = "Arial", height = 3, width = 4)
# my_comparisons1 <- list(c("Up-Up", "Down-Down"), c("Up-Up", "Up-Down"), c("Down-Down", "Up-Down"))
par(xpd=NA)
my_comparisons2 <- list(c("CCs", "HHs"), c("HHs", "CHs"), c("CCs", "CHs"))
ggplot(chrdis_df, aes(x=Region_Type, y = chr_dist, fill = Region_Type)) + 
    geom_violin(adjust = .5) +
    geom_boxplot(width = .1, fill = "black", outlier.color = NA) + 
    stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) +
    stat_compare_means(comparisons = my_comparisons2, label = "p.signif", method = "wilcox.test", size = 3)+
    scale_fill_manual(values = mycolor) + 
    scale_color_manual(values = mycolor) +
    labs(x="",y="Distance", fill = "") + 
    theme(axis.title=element_text(size=8,face="bold"),
          axis.text=element_text(size=8),
          axis.text.x = element_text(face = "bold"),     # adjust the x label with 45 angle
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(face="bold"),
          legend.position = "none",                      # no legend
          panel.grid.major  = element_blank(),           # remove backgroud 
          panel.background = element_blank(),            # remove backgroud 
          axis.line=element_line(size = .3, colour="black"), 
          strip.text.x=element_text(size=8,face="bold"))
          # plot.margin = unit(c(20, 3, 3, 3), "points"))
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distance_tissuesplit.pdf")

cc <- chrdis_df$chr_dist[chrdis_df$Region_Type=="CCs"]
hh <- chrdis_df$chr_dist[chrdis_df$Region_Type=="HHs"]
ch <- chrdis_df$chr_dist[chrdis_df$Region_Type=="CHs"]

tmp <- data.frame(Median = sapply(split(chrdis_df, f = chrdis_df$Region_Type), function(x) {median(x$chr_dist)}), 
                  Mean = sapply(split(chrdis_df, f = chrdis_df$Region_Type), function(x) {mean(x$chr_dist)}))
write.table(tmp, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distance_tissuesplit.txt", quote = F, sep = "\t")

tmp <- data.frame(cc.hh = wilcox.test(cc, hh)$p.value, 
                  hh.ch = wilcox.test(hh, ch)$p.value, 
                  cc.ch = wilcox.test(cc, ch)$p.value)
write.table(tmp, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distance_p_tissuesplit.txt", quote = F, sep = "\t", row.names = F)

# ----------------------------------------------------------------------------------------------------
# V. mini volcano plot for fig1
# ----------------------------------------------------------------------------------------------------
tmp <- load("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/RData/Brain_Caudate.RData")
predata <- data.frame(DESeq2::results(tdds, alpha = 0.05))
library(ggplot2)

plotDF <- predata[,c("log2FoldChange","padj")]
plotDF$regulation <- "no"
plotDF$regulation[plotDF$padj < 0.05 & plotDF$log2FoldChange >0] <- "up"
plotDF$regulation[plotDF$padj < 0.05 & plotDF$log2FoldChange <0] <- "down"

downColor <- "red"
noChangeColor <- "black"
upColor <- "green4"
plotDF <- 




xlimit <- round(max(abs(plotDF$log2FoldChange)))
p <- ggplot(plotDF, aes(x = log2FoldChange, y = -1 * log10(padj), col=regulation))
p <- p + geom_point(size=0.8)
p <- p + scale_color_manual(name = "Regulation", values = c(downColor,noChangeColor,upColor))
p <- p + scale_shape_manual(name = "Regulation" , values = c(20,20,20))
p <- p + scale_x_continuous(limits = c(-2, 2.5))
#p <- p + labs(subtitle = plotSubTitle,x = "log2(FoldChange)", y = "-log10(Adjust pvalue)")
p <- p + theme(legend.text = element_text(size = 10), 
               legend.title = element_text(size = 11), 
               axis.title = element_text(size = 10), 
               panel.grid.major  = element_blank(),           # remove backgroud 
               panel.background = element_blank(),            # remove backgroud 
               axis.line=element_line(size = .3, colour="black"))
pdf("/data/MyProgram/Final_diabrain/4.plots/WGCNA_Caudate/mini_volcano.pdf", width = 12, height = 8)
print(p)
dev.off()
