# ====================================================================================================
# Author:   Aaron
# Function:  plot the diabetes vs ctrl sample counts in tissues 
# Version:  1.0
# Date:     Feb 26, 2019
# ====================================================================================================

library(ggplot2)
library(extrafont)
rm(list=ls())
loadfonts()
options(stringsAsFactors = FALSE)
setwd("/data/MyProgram/Final_diabrain/4.plots/Sample_size/")
sam <- read.table("/data/MyProgram/Final_diabrain/2.results/Cov_optimal/sample_size.tab", header = TRUE)
sam$Tissue <- sub("Brain_", "", sam$Tissue)
pre_data <- rbind(data.frame(Count = sam$Ctrl, Tissue = sam$Tissue, Type = rep("Ctrl", 13)), 
                  data.frame(Count = sam$T2D, Tissue = sam$Tissue, Type = rep("T2D", 13)))

pre_data$Tissue <- gsub("_", " ", pre_data$Tissue)

pdf("./sample_distribute.pdf", width = 8, height = 6)
ggplot(pre_data, aes(x=Type, y = Count, group = Type, color = Type, fill = Type)) + 
    geom_bar(stat = "identity", width = 0.7) + 
    geom_text(size = 3, aes(x=Type, y=Count+5, label = Count, color = Type)) + # add number of each bin
    facet_wrap(~Tissue,ncol=3) +
    coord_flip() + 
    labs(x="Sample Type",y="Count")+
    theme(axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_text(size=8,face="bold"),
          legend.title=element_text(size=8,face="bold"),
          legend.text=element_text(size=8, face="bold"),
          strip.text.x=element_text(size=8,face="bold"))
dev.off()
# embed_fonts("/data/MyProgram/Final_diabrain/4.plots/Sample_size/sample_distribute.pdf")
