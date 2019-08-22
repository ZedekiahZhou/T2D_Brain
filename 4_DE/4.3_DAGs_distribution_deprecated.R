# ====================================================================================================
# Author:   Aaron
# Function: DEGs distribution
# Version:  1.0
# Date:     Mar 18, 2019
# ====================================================================================================
library(ggplot2)
library(ggpubr)
library(scales)
rm(list = ls())
inputfiles <- dir("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/result05/res/", full.names = T)
options(stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------------------------------
# I. get the DEGs from all the 13 tissue
# ----------------------------------------------------------------------------------------------------
upgene <- character() -> downgene
for (i in 1:length(inputfiles)) {
    tmp <- read.table(inputfiles[i], header = T)
    upgene <- c(upgene, tmp$Gene[tmp$log2FoldChange > 0])
    downgene <- c(downgene, tmp$Gene[tmp$log2FoldChange < 0])
}

upgene <- unique(upgene)
downgene <- unique(downgene)

# get the annotations of all genes
gene_info <- read.table("/data/MyProgram/Final_diabrain/1.clean/gen_anno/Total_Gene_Info2.list", header = T)
upgene_info <- gene_info[match(upgene, gene_info$ID), ]
downgene_info <- gene_info[match(downgene, gene_info$ID), ]

# ----------------------------------------------------------------------------------------------------
# II. Gene type
# ----------------------------------------------------------------------------------------------------

# prepare for bar plot
pre_data <- rbind(data.frame(upgene_info, Regulation = "Up"), data.frame(downgene_info, Regulation = "Down"))
pre_data$GeneType <- as.factor(pre_data$GeneType)

# summary for the genetype counts
upgene_type <- as.data.frame(table(pre_data$GeneType[pre_data$Regulation == "Up"]))
downgene_type <- as.data.frame(table(pre_data$GeneType[pre_data$Regulation == "Down"]))
colnames(upgene_type) <- c("GeneType", "Freq.Up")
colnames(downgene_type) <- c("GeneType", "Freq.Down")
upgene_type$Percent.Up <- percent(upgene_type$Freq.Up/sum(upgene_type$Freq.Up))
downgene_type$Percent.Down <- percent(downgene_type$Freq.Down/sum(downgene_type$Freq.Down))
total_type <- merge(upgene_type, downgene_type, by = "GeneType")
write.table(total_type, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Gene_type_distribution.txt", sep = "\t", quote = F, row.names = F)

genetype_description <- data.frame(Genetype = c("protein_coding", "lincRNA", "antisense", "pseudogene"), 
                                   Description = c("Protein coding", "LincRNA", "Antisense RNA", "Pseudogene"))
idx <- match(pre_data$GeneType, genetype_description$Genetype)
pre_data$GeneType <- ifelse(is.na(idx), "Others", genetype_description$Description[idx])

# bar plot
pdf("/data/MyProgram/Final_diabrain/4.plots/DEGs/Gene_type_distribution.pdf", height = 4, width = 8)
ggplot(pre_data, aes(x=GeneType, y = ..prop.., group = Regulation, color = Regulation, fill = Regulation)) + 
    geom_bar(width = 0.7, position = "dodge") +      # plot groups side by side
    coord_flip() +                 # reverse x and y axis
    scale_y_continuous(labels = percent, expand = c(0, 0)) +    # remove the gab between geom and x axis
    labs(x="Gene Type",y="Proportion")+
    theme(axis.title.x=element_text(size=15,face="bold"),
          axis.title.y=element_text(size=15,face="bold"),
          axis.text=element_text(size=12),
          axis.text.y = element_text(face = "bold"),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(face="bold"),
          panel.grid.major  = element_blank(),            # remove backgroud 
          panel.background = element_blank(),             # remove backgroud 
          axis.line=element_line(size = .1, colour="black"),  # add axis
          strip.text.x=element_text(size=10,face="bold"))
dev.off()



# ----------------------------------------------------------------------------------------------------
# III. chromosome distribution
# ----------------------------------------------------------------------------------------------------
pre_data$CHR <- factor(pre_data$CHR, levels = c(as.character(1:22), "X", "Y", "MT"))

# summary for the chromosome counts
upgene_chr <- as.data.frame(table(pre_data$CHR[pre_data$Regulation == "Up"]))
downgene_chr <- as.data.frame(table(pre_data$CHR[pre_data$Regulation == "Down"]))
colnames(upgene_chr) <- c("Chr", "Freq.Up")
colnames(downgene_chr) <- c("Chr", "Freq.Down")
upgene_chr$Percent.Up <- percent(upgene_chr$Freq.Up/sum(upgene_chr$Freq.Up))
downgene_chr$Percent.Down <- percent(downgene_chr$Freq.Down/sum(downgene_chr$Freq.Down))
total_chr <- merge(upgene_chr, downgene_chr, by = "Chr")
write.table(total_chr, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distribution.txt", sep = "\t", quote = F, row.names = F)

pdf("/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distribution.pdf", height = 4, width = 16)
ggplot(pre_data, aes(x=CHR, y = ..prop.., group = Regulation, color = Regulation, fill = Regulation)) + 
    geom_bar(width = 0.5, position = "dodge") + 
    # coord_flip() + 
    labs(x="Chromosome",y="Proportion")+
    scale_y_continuous(labels = percent, expand = c(0, 0)) + 
    theme(axis.title.x=element_text(size=15,face="bold"),
          axis.title.y=element_text(size=15,face="bold"),
          axis.text=element_text(size=12),
          axis.text.x = element_text(face = "bold"),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(face="bold"),
          panel.grid.major  = element_blank(),
          panel.background = element_blank(), 
          axis.line=element_line(size = .1, colour="black"), 
          strip.text.x=element_text(size=10,face="bold"))
dev.off()

# ----------------------------------------------------------------------------------------------------
# IV. chromosome distance
# ----------------------------------------------------------------------------------------------------
# pre_data <- pre_data[order(pre_data$CHR), ]
chr_split <- split(pre_data, f = pre_data$CHR)
chr_dist = as.integer()
Gene_Pair = as.character()
Pair_Type = as.character()
l <- 0

for (i in 1:length(chr_split)) {
    if (dim(chr_split[[i]])[1] == 1) next
    for (j in 1:(dim(chr_split[[i]])[1]-1) ) {
        for (k in (j+1):dim(chr_split[[i]])[1] ) {
            dis1 <- chr_split[[i]]$Start[j] - chr_split[[i]]$End[k]
            dis2 <- chr_split[[i]]$Start[k] - chr_split[[i]]$End[j]
            l <- l+1
            chr_dist[l] <- max(dis1, dis2, 0)
            Gene_Pair[l] <- paste(chr_split[[i]]$ID[j], chr_split[[i]]$ID[k], sep = "-")
            Pair_Type[l] <- paste(chr_split[[i]]$Regulation[j], chr_split[[i]]$Regulation[k], sep = "-")
        }
    }
}

chrdis_df <- data.frame(chr_dist, Gene_Pair, Pair_Type)

pdf("/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distance.pdf", height = 6, width = 8)
my_comparisons <- list(c("Up-Up", "Down-Down"), c("Up-Up", "Up-Down"), c("Down-Down", "Up-Down"))
ggplot(chrdis_df, aes(x=Pair_Type, y = chr_dist, fill = Pair_Type)) + 
    geom_violin(adjust = .5) +
    geom_boxplot(width = .1, fill = "black", outlier.color = NA) + 
    stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")+
    # coord_flip() + 
    labs(x="",y="Distance", fill = "") + 
    theme(axis.title.x=element_text(size=15,face="bold"),
          axis.title.y=element_text(size=15,face="bold"),
          axis.text=element_text(size=12),
          axis.text.x = element_text(face = "bold"),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(face="bold"),
          legend.position = "none", 
          panel.grid.major  = element_blank(),
          panel.background = element_blank(), 
          axis.line=element_line(size = .1, colour="black"), 
          strip.text.x=element_text(size=10,face="bold"))
dev.off()

uu <- chrdis_df$chr_dist[chrdis_df$Pair_Type=="Up-Up"]
dd <- chrdis_df$chr_dist[chrdis_df$Pair_Type=="Down-Down"]
ud <- chrdis_df$chr_dist[chrdis_df$Pair_Type=="Up-Down"]

tmp <- data.frame(Median = sapply(split(chrdis_df, f = chrdis_df$Pair_Type), function(x) {median(x$chr_dist)}), 
                     Mean = sapply(split(chrdis_df, f = chrdis_df$Pair_Type), function(x) {mean(x$chr_dist)}))
write.table(tmp, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distance.txt", quote = F, sep = "\t")

tmp <- data.frame(uu.dd = wilcox.test(uu, dd)$p.value, 
                  uu.ud = wilcox.test(uu, ud)$p.value, 
                  ud.dd = wilcox.test(ud, dd)$p.value)
write.table(tmp, "/data/MyProgram/Final_diabrain/4.plots/DEGs/Chromosome_distance_p.txt", quote = F, sep = "\t", row.names = F)
